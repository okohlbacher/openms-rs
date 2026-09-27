// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PTMXML files: a flat list of post-translational modifications.
//!
//! Ports `FORMAT/PTMXMLFile.h` and its handler `FORMAT/HANDLERS/PTMXMLHandler.h`.
//! A document is
//!
//! ```text
//! <PTMs>
//!     <PTM>
//!         <name>TEST</name>
//!         <composition>N2O2-CH3</composition>
//!         <possible_amino_acids>KLR</possible_amino_acids>
//!     </PTM>
//! </PTMs>
//! ```
//!
//! and loads into a map from the modification name to its composition and
//! the one-letter codes of the residues it may modify, [`PtmInformation`].
//! See `docs/PTM_XML_SUPPORT.md` for the member mapping and the evidence.

use super::sax_xml::{SaxEvent, SaxLimits, parse};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::Path;

/// PTM information keyed by modification name: the pair is
/// `(composition, possible_amino_acids)`, both exactly as the file holds them
/// after trimming.
///
/// Source `std::map<std::string, std::pair<std::string, std::string>>`; a
/// [`BTreeMap`] keeps the byte order `std::map` iterates in, which is also the
/// order [`write`] emits.
pub type PtmInformation = BTreeMap<String, (String, String)>;

/// Resource bounds for reading one PTMXML document. Exceeding any of them is
/// [`Error::InvalidValue`] and yields no result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Decompressed input bytes.
    pub max_input_bytes: u64,
    /// Open elements at any one time.
    pub max_depth: usize,
    /// Parser events in the document.
    pub max_events: u64,
    /// Bytes of one run of character data.
    pub max_text_bytes: usize,
    /// Distinct modification names in the result.
    pub max_entries: usize,
    /// Bytes of all names, compositions and residue lists in the result.
    ///
    /// The source carries the last `<composition>` over to every later
    /// modification that lacks one, so a small document can otherwise ask
    /// for a large result; this bounds that amplification.
    pub max_output_bytes: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_depth: 64,
            max_events: 50_000_000,
            max_text_bytes: 16 * 1024 * 1024,
            max_entries: 1_000_000,
            max_output_bytes: 512 * 1024 * 1024,
        }
    }
}

impl ReadLimits {
    fn sax(&self) -> SaxLimits {
        SaxLimits {
            max_input_bytes: self.max_input_bytes,
            max_depth: self.max_depth,
            max_events: self.max_events,
            max_attributes: 256,
            max_text_bytes: self.max_text_bytes,
        }
    }
}

/// Source `StringUtils::trimmed`: space, tab, line feed and carriage return.
fn trimmed(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\n', '\r'])
}

/// The state of source `PTMXMLHandler`, member for member.
#[derive(Default)]
struct Handler {
    /// `tag_`: the trimmed name of the element most recently opened, cleared
    /// by any end tag.
    tag: String,
    /// `name_`: the last `<name>` text; it survives into later `<PTM>`s.
    name: String,
    /// `composition_`: the last `<composition>` text; it too survives.
    composition: String,
    /// `open_tag_`: an element has opened and no element has closed since.
    open: bool,
    output: PtmInformation,
    output_bytes: usize,
}

impl Handler {
    fn event(&mut self, event: SaxEvent<'_>, limits: &ReadLimits) -> Result<()> {
        match event {
            SaxEvent::Start(qname, _) => {
                self.tag = trimmed(qname).to_owned();
                self.open = true;
            }
            SaxEvent::End(_) => {
                self.tag.clear();
                self.open = false;
            }
            SaxEvent::Characters(chars) => {
                if !self.open {
                    return Ok(());
                }
                match self.tag.as_str() {
                    "name" => self.name = trimmed(chars).to_owned(),
                    "composition" => self.composition = trimmed(chars).to_owned(),
                    "possible_amino_acids" => {
                        self.insert(trimmed(chars).to_owned(), limits)?;
                    }
                    _ => {}
                }
            }
        }
        Ok(())
    }

    /// `ptm_informations_[name_] = make_pair(composition_, amino_acids)`,
    /// with the result bounded before it grows.
    fn insert(&mut self, amino_acids: String, limits: &ReadLimits) -> Result<()> {
        let added = self
            .name
            .len()
            .saturating_add(self.composition.len())
            .saturating_add(amino_acids.len());
        let removed = self
            .output
            .get(&self.name)
            .map_or(0, |(c, a)| self.name.len() + c.len() + a.len());
        let bytes = self.output_bytes.saturating_sub(removed).saturating_add(added);
        if bytes > limits.max_output_bytes {
            return Err(Error::InvalidValue("PTMXML output byte limit exceeded".into()));
        }
        if !self.output.contains_key(&self.name) && self.output.len() == limits.max_entries {
            return Err(Error::InvalidValue("PTMXML entry limit exceeded".into()));
        }
        self.output
            .insert(self.name.clone(), (self.composition.clone(), amino_acids));
        self.output_bytes = bytes;
        Ok(())
    }
}

/// Read a PTMXML document with the default [`ReadLimits`].
///
/// Source `PTMXMLFile::load`, which clears its output map first; this returns
/// a new map instead, and only on success.
///
/// The reader is the source handler's state machine, so a document that is
/// not shaped like the example above is read the way the source reads it:
///
/// - Only character data directly inside `<name>`, `<composition>` and
///   `<possible_amino_acids>` matters, and element nesting is not checked:
///   the three may appear under any parent, in any order.
/// - An entry is stored when `<possible_amino_acids>` has text, under the
///   most recent `<name>` and with the most recent `<composition>`, both
///   carried over from earlier `<PTM>`s when the current one lacks them. An
///   empty `<possible_amino_acids/>` stores nothing, and a later entry with
///   the same name replaces an earlier one.
/// - Every text is trimmed of spaces, tabs and line ends.
/// - Character data interrupted by a CDATA section, comment or processing
///   instruction arrives in pieces, as Xerces delivers it, and each piece
///   overwrites the last.
///
/// # Errors
///
/// [`Error::Parse`] when the document is not well-formed XML, which is the
/// source's `Exception::ParseError`; [`Error::Unsupported`] for a document
/// type declaration or a non-UTF-8 encoding; [`Error::InvalidValue`] when a
/// [`ReadLimits`] bound is exceeded; [`Error::Io`] when the input fails.
pub fn read(input: impl BufRead) -> Result<PtmInformation> {
    read_with_limits(input, &ReadLimits::default())
}

/// [`read`] under explicit resource bounds.
///
/// # Errors
///
/// As [`read`].
pub fn read_with_limits(input: impl BufRead, limits: &ReadLimits) -> Result<PtmInformation> {
    let mut handler = Handler::default();
    parse(input, &limits.sax(), |event| handler.event(event, limits))?;
    Ok(handler.output)
}

/// Load a PTMXML file with the default [`ReadLimits`]; gzip and bzip2 input
/// are recognised by their leading bytes, as the source's `XMLFile::parse_`
/// does.
///
/// # Errors
///
/// [`Error::Io`] when the file cannot be opened, which is the source's
/// `Exception::FileNotFound`; otherwise as [`read`].
pub fn load(path: impl AsRef<Path>) -> Result<PtmInformation> {
    load_with_limits(path, &ReadLimits::default())
}

/// [`load`] under explicit resource bounds.
///
/// # Errors
///
/// As [`load`].
pub fn load_with_limits(path: impl AsRef<Path>, limits: &ReadLimits) -> Result<PtmInformation> {
    read_with_limits(super::path_io::open(path.as_ref())?, limits)
}

/// Escape the five characters the source's `XMLHandler::writeXMLEscape`
/// replaces, refusing characters XML 1.0 cannot carry.
fn escaped(text: &str, output: &mut String) -> Result<()> {
    for c in text.chars() {
        match c {
            '&' => output.push_str("&amp;"),
            '<' => output.push_str("&lt;"),
            '>' => output.push_str("&gt;"),
            '"' => output.push_str("&quot;"),
            '\'' => output.push_str("&apos;"),
            '\t' | '\n' | '\r' => output.push(c),
            c if c.is_control() || matches!(c, '\u{FFFE}' | '\u{FFFF}') => {
                return Err(Error::InvalidValue(format!(
                    "PTMXML cannot represent the character U+{:04X}",
                    u32::from(c)
                )));
            }
            c => output.push(c),
        }
    }
    Ok(())
}

/// Serialize `ptm_informations` exactly as source `PTMXMLHandler::writeTo`
/// does: no XML declaration, tab indentation, `\n` line ends, one `<PTM>`
/// per entry in key order.
///
/// Source `writeTo` inserts the strings unescaped, so a name holding `&` or
/// `<` produces a file no XML parser reads back. This escapes the five
/// characters `XMLHandler::writeXMLEscape` escapes; for any other text the
/// bytes are the source's.
///
/// The document is built completely before anything is written.
///
/// # Errors
///
/// [`Error::InvalidValue`] when a string holds a control character other than
/// tab, line feed or carriage return, which XML 1.0 cannot carry;
/// [`Error::Io`] when the output fails.
pub fn write(mut output: impl Write, ptm_informations: &PtmInformation) -> Result<()> {
    let text = document(ptm_informations)?;
    output.write_all(text.as_bytes())?;
    output.flush()?;
    Ok(())
}

fn document(ptm_informations: &PtmInformation) -> Result<String> {
    let mut text = String::from("<PTMs>\n");
    for (name, (composition, amino_acids)) in ptm_informations {
        text.push_str("\t<PTM>\n\t\t<name>");
        escaped(name, &mut text)?;
        text.push_str("</name>\n\t\t<composition>");
        escaped(composition, &mut text)?;
        text.push_str("</composition>\n\t\t<possible_amino_acids>");
        escaped(amino_acids, &mut text)?;
        text.push_str("</possible_amino_acids>\n\t</PTM>\n");
    }
    text.push_str("</PTMs>\n");
    Ok(text)
}

/// Store `ptm_informations` as a PTMXML file.
///
/// Source `PTMXMLFile::store` through `XMLFile::save_`: a name ending in
/// `.gz` or `.bz2` is written compressed. The file is published atomically,
/// from a sibling temporary file, only after the whole document has been
/// written; the source truncates the destination first.
///
/// # Errors
///
/// As [`write`]; [`Error::Io`] when the file cannot be created, which is the
/// source's `Exception::UnableToCreateFile`.
pub fn store(path: impl AsRef<Path>, ptm_informations: &PtmInformation) -> Result<()> {
    let text = document(ptm_informations)?;
    super::path_io::write(path.as_ref(), |writer| {
        writer.write_all(text.as_bytes())?;
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn source_state_carries_over_between_entries() {
        let text = "<PTMs><PTM><name>A</name><composition>C1</composition>\
                    <possible_amino_acids> K </possible_amino_acids></PTM>\
                    <PTM><name>B</name><possible_amino_acids>R</possible_amino_acids></PTM>\
                    <PTM><name>C</name><possible_amino_acids/></PTM></PTMs>";
        let map = read(text.as_bytes()).unwrap();
        assert_eq!(map.len(), 2);
        assert_eq!(map["A"], ("C1".into(), "K".into()));
        // B has no composition of its own and inherits A's, as in the source.
        assert_eq!(map["B"], ("C1".into(), "R".into()));
    }

    #[test]
    fn output_bound_stops_carried_composition_amplification() {
        let limits = ReadLimits {
            max_output_bytes: 8,
            ..ReadLimits::default()
        };
        let text = "<PTMs><name>A</name><composition>CCCC</composition>\
                    <possible_amino_acids>K</possible_amino_acids>\
                    <name>B</name><possible_amino_acids>K</possible_amino_acids></PTMs>";
        assert!(matches!(
            read_with_limits(text.as_bytes(), &limits),
            Err(Error::InvalidValue(_))
        ));
        let replaced = "<PTMs><name>A</name><composition>CCCC</composition>\
                        <possible_amino_acids>K</possible_amino_acids>\
                        <possible_amino_acids>R</possible_amino_acids></PTMs>";
        let map = read_with_limits(replaced.as_bytes(), &limits).unwrap();
        assert_eq!(map["A"], ("CCCC".into(), "R".into()));
    }

    #[test]
    fn writer_refuses_characters_xml_cannot_carry() {
        let map = PtmInformation::from([("a\u{1}".to_owned(), (String::new(), "K".to_owned()))]);
        assert!(matches!(
            write(Vec::new(), &map),
            Err(Error::InvalidValue(_))
        ));
    }
}
