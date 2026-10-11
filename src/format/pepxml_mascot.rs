// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for the pepXML that Mascot exports, ported from
//! `FORMAT/PepXMLFileMascot.h`.
//!
//! A schema for this format is Matrix Science's
//! `http://www.matrixscience.com/xmlns/schema/pepXML_v18/pepXML_v18.xsd`; it is
//! not consulted, as the source bundles no XSD for it either.
//!
//! The reader keeps only peptide sequences with their modifications, keyed by
//! the `spectrum` attribute of their `<spectrum_query>`. Scores, proteins and
//! every other element are ignored; for full pepXML identifications use
//! [`crate::format::pepxml`], the port of `PepXMLFile.h`.
//!
//! How a sequence is built, as the source builds it at the end of each
//! `<search_hit>`:
//!
//! 1. The `peptide` attribute is parsed as an [`AASequence`](crate::chemistry::AASequence).
//! 2. Each `<mod_aminoacid_mass>` of the hit is matched by mass to the first
//!    *variable* `<aminoacid_modification>` or `<terminal_modification>` of
//!    the search summary whose `mass` is within [`MODIFICATION_TOLERANCE`](crate::format::pepxml_mascot::MODIFICATION_TOLERANCE)
//!    (strictly less than 0.001 Da). Its description (for example
//!    `Oxidation (M)`) is applied: a description containing `C-term` or
//!    `N-term` is also set, whole, as the C- or N-terminal modification; a
//!    description of exactly two space-separated words sets the first word as
//!    the modification of the residue at the 1-based `position`. Any other
//!    description, and a mass that matches nothing (an empty description), is
//!    the source's non-fatal error `Cannot parse modification '<description>@<position>'`,
//!    logged, and the hit is kept.
//! 3. Every *fixed* modification (any `variable` other than `Y`) of two words is
//!    applied: `(C-term)` and `(N-term)` set the terminal modification, and any
//!    other second word, with its parentheses removed, lists the residues the
//!    first word is set on (`Phospho (STY)` sets `Phospho` on every S, T and
//!    Y). A fixed modification of any other shape is the non-fatal error
//!    `Cannot parse fixed modification '<description>'`, logged.
//!
//! The sequences of the hits of one `<spectrum_query>` are stored, in document
//! order, under its title when the query ends; a title seen before keeps its
//! first entry, as `std::map::insert` does. The API mapping, preserved source
//! conventions and every native difference are listed in
//! `docs/PEPXML_MASCOT_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::data_structures::string_utils::to_double;
use crate::format::sax_xml::{
    self, SaxAttributes as Attributes, SaxEncoding, SaxEvent, SaxLimits, xerces_parse_int,
};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Mass window of the source's `matchModification_`, in Da: a variable
/// modification matches when its `mass` differs from the observed mass by
/// strictly less than this.
pub const MODIFICATION_TOLERANCE: f64 = 0.001;

/// Resource ceilings of one read, each checked before the work it bounds.
///
/// The defaults accept every upstream fixture and any Mascot export of a
/// single search.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum document size in bytes.
    pub max_bytes: usize,
    /// Maximum XML events, including those of ignored elements.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum `<aminoacid_modification>` and `<terminal_modification>`
    /// definitions, fixed and variable together. Every observed modification
    /// mass is compared with every variable definition, so this bounds that
    /// work.
    pub max_modification_definitions: usize,
    /// Maximum `<mod_aminoacid_mass>` elements within one `<search_hit>`.
    pub max_modifications_per_hit: usize,
    /// Maximum residues of one `peptide` attribute. A fixed modification is
    /// set residue by residue, and each setting rebuilds the sequence, so the
    /// work of one hit grows with the square of its length.
    pub max_residues: usize,
    /// Maximum sequences held, over all spectrum titles and the pending
    /// query.
    pub max_sequences: usize,
    /// Maximum retained log lines; later ones are counted in
    /// [`MascotPepXmlDocument::suppressed_log_lines`].
    pub max_log_lines: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1 << 30,
            max_events: 100_000_000,
            max_depth: 64,
            max_modification_definitions: 1_024,
            max_modifications_per_hit: 1_024,
            max_residues: 10_000,
            max_sequences: 10_000_000,
            max_log_lines: 10_000,
        }
    }
}

/// What one read produced: the source's `peptides` out-parameter and its log.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MascotPepXmlDocument {
    /// The peptide sequences of each spectrum, keyed by the `spectrum`
    /// attribute of its `<spectrum_query>`, in byte order of the title (the
    /// order of the source's `std::map<std::string, ...>`); the sequences of
    /// one query are in document order.
    pub peptides: BTreeMap<String, Vec<AASequence>>,
    /// The source's non-fatal errors (`XMLHandler::error`, which it writes to
    /// `OPENMS_LOG_ERROR`), in order, each as
    /// `Non-fatal error while loading '<file>': <message>`.
    pub log: Vec<String>,
    /// Log lines beyond [`ReadLimits::max_log_lines`].
    pub suppressed_log_lines: usize,
}

/// The source's `PepXMLFileMascot`, a stateless reader.
///
/// The source keeps the modification definitions and the hit under
/// construction in members, resets them only after a successful load, and
/// never resets the sequences of a query that has not ended; a load after one
/// that threw therefore still matches masses against the failed file's
/// variable modifications (see `docs/PEPXML_MASCOT_SUPPORT.md`). Every read
/// here starts from empty state.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PepXmlFileMascot;

impl PepXmlFileMascot {
    /// A reader, the source's constructor.
    pub fn new() -> Self {
        Self
    }

    /// Load peptide sequences with modifications out of a Mascot pepXML file,
    /// as `load(filename, peptides)`.
    ///
    /// # Errors
    ///
    /// - [`Error::Io`] when the file cannot be opened or read (the source's
    ///   `Exception::FileNotFound` for a missing file), and
    ///   [`Error::InvalidValue`] when a [`ReadLimits`] ceiling is exceeded
    ///   (except [`ReadLimits::max_bytes`], which the encoding-aware document
    ///   decoder enforces and reports as [`Error::Parse`]);
    /// - [`Error::Parse`] (the source's `Exception::ParseError`) for a
    ///   document that is not well-formed XML, a missing required attribute
    ///   (`variable`, `description` and, for a variable modification, `mass`
    ///   of a modification definition; `spectrum`; `peptide`; `position` and
    ///   `mass` of `<mod_aminoacid_mass>`), a `position` that is not an
    ///   integer, a `mass` that is not a number, or a `peptide` that does not
    ///   parse as a sequence (the source's `AASequence::fromString` error);
    /// - [`Error::InvalidRange`] when a `<mod_aminoacid_mass>` position is
    ///   not a residue of its peptide (the source's
    ///   `Exception::IndexOverflow`; position 0 wraps to the largest index
    ///   there);
    /// - [`Error::InvalidValue`] when a matched or fixed modification names nothing the modification registry has for
    ///   that residue or terminus (the source's
    ///   `Exception::ElementNotFound`);
    /// - [`Error::Unsupported`] for a DTD.
    ///
    /// Stricter than the source: a `<mod_aminoacid_mass>` position outside the
    /// 32-bit range is refused (the source's `(int)` cast wraps it), and a
    /// terminal modification on an empty `peptide` is refused
    /// ([`AASequence`] has no terminus to hold it).
    pub fn load(&self, path: impl AsRef<Path>) -> Result<MascotPepXmlDocument> {
        self.load_with_limits(path, &ReadLimits::default())
    }

    /// [`Self::load`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn load_with_limits(
        &self,
        path: impl AsRef<Path>,
        limits: &ReadLimits,
    ) -> Result<MascotPepXmlDocument> {
        let path = path.as_ref();
        let file = File::open(path)?;
        read_document(BufReader::new(file), &path.display().to_string(), limits)
    }

    /// Read Mascot pepXML from `reader`; see [`Self::load`]. The file name in
    /// log lines is empty.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read(&self, reader: impl BufRead) -> Result<MascotPepXmlDocument> {
        read_document(reader, "", &ReadLimits::default())
    }

    /// [`Self::read`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read_with_limits(
        &self,
        reader: impl BufRead,
        limits: &ReadLimits,
    ) -> Result<MascotPepXmlDocument> {
        read_document(reader, "", limits)
    }
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn limit(what: &str) -> Error {
    Error::InvalidValue(format!("Mascot pepXML {what} limit exceeded"))
}

/// `attributeAsDouble_`: `StringUtils::toDouble`, whose conversion error ends
/// the parse.
fn double(attributes: &Attributes, name: &str) -> Result<f64> {
    let text = attributes.required(name)?;
    to_double(text).map_err(|_| {
        parse_error(format!(
            "Could not convert '{text}' in '{name}' to a double"
        ))
    })
}

fn read_document(
    reader: impl BufRead,
    file: &str,
    limits: &ReadLimits,
) -> Result<MascotPepXmlDocument> {
    let sax = SaxLimits {
        max_input_bytes: u64::try_from(limits.max_bytes).unwrap_or(u64::MAX),
        max_events: u64::try_from(limits.max_events).unwrap_or(u64::MAX),
        max_depth: limits.max_depth,
        max_attributes: 256,
        // Character data is not read; a run cannot exceed the document.
        max_text_bytes: limits.max_bytes.max(1),
    };
    let mut handler = Handler {
        file,
        limits: *limits,
        document: MascotPepXmlDocument::default(),
        title: String::new(),
        sequence: String::new(),
        modifications: Vec::new(),
        pending: Vec::new(),
        sequences: 0,
        fixed: Vec::new(),
        variable: Vec::new(),
    };
    sax_xml::parse_with(reader, SaxEncoding::Declared, &sax, |event| match event {
        SaxEvent::Start(name, attributes) => handler.start(name, attributes),
        SaxEvent::End(name) => handler.end(name),
        SaxEvent::Characters(_) => Ok(()),
    })?;
    Ok(handler.document)
}

/// The handler state of the source's `PepXMLFileMascot`: `actual_title_`,
/// `actual_sequence_`, `actual_modifications_`, `actual_aa_sequences_`,
/// `fixed_modifications_` and `variable_modifications_`, with `peptides_` as
/// the document.
struct Handler<'a> {
    file: &'a str,
    limits: ReadLimits,
    document: MascotPepXmlDocument,
    /// `actual_title_`: the `spectrum` of the last `<spectrum_query>` opened.
    title: String,
    /// `actual_sequence_`: the `peptide` of the last `<search_hit>` opened.
    sequence: String,
    /// `actual_modifications_`: description (empty when no mass matched) and
    /// the 1-based position, as the source's `UInt`.
    modifications: Vec<(String, u32)>,
    /// `actual_aa_sequences_`: the hits of the query not yet ended. A hit
    /// outside any query stays here and is stored with the next query, as in
    /// the source.
    pending: Vec<AASequence>,
    /// Sequences held in the document and `pending`.
    sequences: usize,
    /// `fixed_modifications_`: descriptions.
    fixed: Vec<String>,
    /// `variable_modifications_`: description and mass.
    variable: Vec<(String, f64)>,
}

impl Handler<'_> {
    /// `XMLHandler::error(LOAD, message)`: logged, not thrown.
    fn error(&mut self, message: String) {
        if self.document.log.len() < self.limits.max_log_lines {
            let line = format!("Non-fatal error while loading '{}': {message}", self.file);
            self.document.log.push(line);
        } else {
            self.document.suppressed_log_lines += 1;
        }
    }

    /// `matchModification_`: the description of the first variable
    /// modification within [`MODIFICATION_TOLERANCE`], or empty.
    fn match_modification(&self, mass: f64) -> String {
        self.variable
            .iter()
            .find(|(_, candidate)| (candidate - mass).abs() < MODIFICATION_TOLERANCE)
            .map(|(description, _)| description.clone())
            .unwrap_or_default()
    }

    fn start(&mut self, tag: &str, attributes: &Attributes) -> Result<()> {
        // Search parameters. The source tests the two definition tags with
        // separate `if`s, the second heading the `else if` chain below.
        if tag == "aminoacid_modification" || tag == "terminal_modification" {
            if self.fixed.len() + self.variable.len() >= self.limits.max_modification_definitions {
                return Err(limit("modification definition"));
            }
            if attributes.required("variable")? == "Y" {
                // The source reads both inside one `emplace_back` call, whose
                // arguments the Release build (GCC) evaluates right to left:
                // with both missing, the error names `mass`.
                let mass = double(attributes, "mass")?;
                let description = attributes.required("description")?.to_owned();
                self.variable.push((description, mass));
            } else {
                self.fixed
                    .push(attributes.required("description")?.to_owned());
            }
            return Ok(());
        }
        match tag {
            "spectrum_query" => {
                self.title = attributes.required("spectrum")?.to_owned();
            }
            "search_hit" => {
                self.sequence = attributes.required("peptide")?.to_owned();
            }
            "mod_aminoacid_mass" => {
                if self.modifications.len() >= self.limits.max_modifications_per_hit {
                    return Err(limit("modifications per hit"));
                }
                // `UInt modification_position = attributeAsInt_(...)`: the
                // signed integer converts to unsigned, so -1 is 4294967295.
                let position = xerces_parse_int(attributes.required("position")?)? as u32;
                let mass = double(attributes, "mass")?;
                let description = self.match_modification(mass);
                self.modifications.push((description, position));
            }
            _ => {}
        }
        Ok(())
    }

    fn end(&mut self, tag: &str) -> Result<()> {
        match tag {
            "search_hit" => {
                if self.sequences >= self.limits.max_sequences {
                    return Err(limit("sequence"));
                }
                let sequence = self.build_hit()?;
                self.pending.push(sequence);
                self.sequences += 1;
                self.modifications.clear();
            }
            "spectrum_query" => {
                let hits = std::mem::take(&mut self.pending);
                if self.document.peptides.contains_key(&self.title) {
                    // `std::map::insert` keeps the first entry for a title.
                    self.sequences -= hits.len();
                } else {
                    self.document.peptides.insert(self.title.clone(), hits);
                }
            }
            _ => {}
        }
        Ok(())
    }

    /// The source's `onEndElement("search_hit")` body.
    fn build_hit(&mut self) -> Result<AASequence> {
        if self.sequence.len() > self.limits.max_residues {
            return Err(limit("peptide residue"));
        }
        // `AASequence::fromString` throws `Exception::ParseError`.
        let mut sequence = AASequence::parse(&permissive(&self.sequence)).map_err(|error| {
            parse_error(format!(
                "Cannot convert string to amino acid sequence: {error} in: {}",
                self.sequence
            ))
        })?;
        let modifications = std::mem::take(&mut self.modifications);
        for (description, position) in &modifications {
            if description.contains("C-term") {
                sequence.set_c_terminal_modification(description)?;
            } else if description.contains("N-term") {
                sequence.set_n_terminal_modification(description)?;
            }
            match two_words(description) {
                Some((name, _)) => {
                    // `it->second - 1` in `UInt`: position 0 wraps.
                    let index = position.wrapping_sub(1) as usize;
                    set_residue(&mut sequence, index, name)?;
                }
                None => {
                    self.error(format!(
                        "Cannot parse modification '{description}@{position}'"
                    ));
                }
            }
        }
        self.modifications = modifications;
        let fixed = std::mem::take(&mut self.fixed);
        let applied = apply_fixed(&mut sequence, &fixed, |message| self.error(message));
        self.fixed = fixed;
        applied?;
        Ok(sequence)
    }
}

fn apply_fixed(
    sequence: &mut AASequence,
    fixed: &[String],
    mut error: impl FnMut(String),
) -> Result<()> {
    for description in fixed {
        let Some((name, origin)) = two_words(description) else {
            error(format!("Cannot parse fixed modification '{description}'"));
            continue;
        };
        if origin == "(C-term)" {
            sequence.set_c_terminal_modification(name)?;
        } else if origin == "(N-term)" {
            sequence.set_n_terminal_modification(name)?;
        } else {
            let origin: String = origin.chars().filter(|c| !matches!(c, '(' | ')')).collect();
            // `hasSubstring(origin, residue.getOneLetterCode())`: the residue
            // letters are ASCII.
            let targets: Vec<usize> = sequence
                .as_str()
                .char_indices()
                .filter(|(_, residue)| origin.contains(*residue))
                .map(|(index, _)| index)
                .collect();
            for index in targets {
                set_residue(sequence, index, name)?;
            }
        }
    }
    Ok(())
}

/// The permissive default of the source's `AASequence::fromString`, which
/// [`AASequence::parse`] does not have: outside round and square brackets a
/// stop codon (`*`, `#`, `+`) is the residue `X` and a space is skipped.
/// Inside brackets, and when the brackets do not balance, the text is left
/// for the parser to accept or refuse.
pub(crate) fn permissive(peptide: &str) -> String {
    let mut depth = 0usize;
    let mut out = String::with_capacity(peptide.len());
    for c in peptide.chars() {
        match c {
            '(' | '[' => {
                depth += 1;
                out.push(c);
            }
            ')' | ']' => {
                if depth == 0 {
                    return peptide.to_owned();
                }
                depth -= 1;
                out.push(c);
            }
            '*' | '#' | '+' if depth == 0 => out.push('X'),
            ' ' if depth == 0 => {}
            _ => out.push(c),
        }
    }
    if depth == 0 { out } else { peptide.to_owned() }
}

/// `StringUtils::split(text, ' ')` yielding exactly two parts; an empty text
/// yields none in the source, so it is not two words either.
fn two_words(text: &str) -> Option<(&str, &str)> {
    let mut parts = text.split(' ');
    match (parts.next(), parts.next(), parts.next()) {
        (Some(first), Some(second), None) if !text.is_empty() => Some((first, second)),
        _ => None,
    }
}

/// `AASequence::setModification(index, name)`, whose index check is the
/// source's `Exception::IndexOverflow`.
fn set_residue(sequence: &mut AASequence, index: usize, name: &str) -> Result<()> {
    if index >= sequence.len() {
        return Err(Error::InvalidRange(format!(
            "Mascot pepXML modification index {index} is outside the {}-residue peptide",
            sequence.len()
        )));
    }
    sequence.set_modification(index, name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_words_follows_the_source_split() {
        assert_eq!(two_words("Oxidation (M)"), Some(("Oxidation", "(M)")));
        assert_eq!(two_words(""), None);
        assert_eq!(two_words("Oxidation"), None);
        assert_eq!(two_words("Acetyl (Protein N-term)"), None);
        // Consecutive separators yield an empty part, as in the source.
        assert_eq!(two_words("a  b"), None);
        assert_eq!(two_words(" x"), Some(("", "x")));
    }
}
