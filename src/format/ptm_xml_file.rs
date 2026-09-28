// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! PTMXML files: the post-translational modification tables the Inspect and
//! Sequest adapters look modification names up in (`Inspect_PTMs.xml`,
//! `Sequest_PTMs.xml`).
//!
//! Ports `FORMAT/PTMXMLFile.h`/`.cpp` and its handler
//! `FORMAT/HANDLERS/PTMXMLHandler.h`/`.cpp`. A document lists one `PTM`
//! element per modification:
//!
//! ```text
//! <PTMs>
//!   <PTM>
//!     <name>TEST</name>
//!     <composition>N2O2-CH3</composition>
//!     <possible_amino_acids>KLR</possible_amino_acids>
//!   </PTM>
//! </PTMs>
//! ```
//!
//! The table maps each name to its composition and its possible residues.

use crate::format::search_engine_xml::{
    drive, read_document, Attributes, SaxHandler, TextEncoding, XmlReadLimits,
};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::Path;

/// A PTM table: modification name to (composition, possible amino acids).
pub type PtmInformations = BTreeMap<String, (String, String)>;

/// PTMXML reader and writer.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PTMXMLFile {
    /// Resource ceilings for [`PTMXMLFile::load`].
    pub limits: XmlReadLimits,
}

impl PTMXMLFile {
    /// A reader with the default resource ceilings.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load the PTM table of a file.
    ///
    /// The source handler's rules are kept: an element's text is trimmed of
    /// space, tab, line feed and carriage return; a `possible_amino_acids`
    /// text stores the entry under the name and composition read last, so
    /// a `PTM` without a `name` or `composition` reuses the previous one, a
    /// repeated name replaces the earlier entry, and element nesting is not
    /// checked.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be opened or read (source
    /// `Exception::FileNotFound`), [`Error::Parse`] when it is not well-formed
    /// XML (source `Exception::ParseError`), [`Error::Unsupported`] for a DTD
    /// or a non-UTF-8 encoding, and [`Error::InvalidValue`] when a limit is
    /// exceeded.
    pub fn load(&self, path: impl AsRef<Path>) -> Result<PtmInformations> {
        self.load_reader(std::fs::File::open(path)?)
    }

    /// Load the PTM table of a document; see [`PTMXMLFile::load`].
    ///
    /// # Errors
    ///
    /// As [`PTMXMLFile::load`].
    pub fn load_reader(&self, reader: impl Read) -> Result<PtmInformations> {
        let text = read_document(reader, TextEncoding::Utf8, &self.limits)?;
        let mut handler = Handler::default();
        drive(&text, TextEncoding::Utf8, &mut handler, &self.limits)?;
        Ok(handler.ptm_informations)
    }

    /// Store a PTM table.
    ///
    /// Writes the source layout (tab indentation, one `PTM` per entry in
    /// name order). Unlike the source, which writes the texts raw, `&`, `<`
    /// and `>` are escaped so that every table round-trips.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be created or written (source
    /// `Exception::UnableToCreateFile`).
    pub fn store(&self, path: impl AsRef<Path>, ptm_informations: &PtmInformations) -> Result<()> {
        let mut file = std::fs::File::create(path)?;
        self.write(&mut file, ptm_informations)?;
        file.flush()?;
        Ok(())
    }

    /// Write a PTM table (source `PTMXMLHandler::writeTo`); see
    /// [`PTMXMLFile::store`].
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when writing fails.
    pub fn write(&self, mut writer: impl Write, ptm_informations: &PtmInformations) -> Result<()> {
        let mut os = String::from("<PTMs>\n");
        for (name, (composition, residues)) in ptm_informations {
            os.push_str("\t<PTM>\n");
            os.push_str(&format!("\t\t<name>{}</name>\n", escape(name)));
            os.push_str(&format!(
                "\t\t<composition>{}</composition>\n",
                escape(composition)
            ));
            os.push_str(&format!(
                "\t\t<possible_amino_acids>{}</possible_amino_acids>\n",
                escape(residues)
            ));
            os.push_str("\t</PTM>\n");
        }
        os.push_str("</PTMs>\n");
        writer.write_all(os.as_bytes()).map_err(Error::from)
    }
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Source `StringUtils::trimmed`: space, tab, line feed, carriage return.
fn trimmed(text: &str) -> String {
    text.trim_matches([' ', '\t', '\n', '\r']).to_owned()
}

#[derive(Default)]
struct Handler {
    ptm_informations: PtmInformations,
    tag: String,
    open_tag: bool,
    name: String,
    composition: String,
}

impl SaxHandler for Handler {
    fn start_element(&mut self, name: &str, _attributes: &Attributes) -> Result<()> {
        self.tag = trimmed(name);
        self.open_tag = true;
        Ok(())
    }

    fn end_element(&mut self, _name: &str) -> Result<()> {
        self.tag.clear();
        self.open_tag = false;
        Ok(())
    }

    fn characters(&mut self, text: &str) -> Result<()> {
        if self.open_tag {
            match self.tag.as_str() {
                "name" => self.name = trimmed(text),
                "composition" => self.composition = trimmed(text),
                "possible_amino_acids" => {
                    self.ptm_informations
                        .insert(self.name.clone(), (self.composition.clone(), trimmed(text)));
                }
                _ => {}
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_table_round_trips_with_escaping() {
        let mut table = PtmInformations::new();
        table.insert("A&B".into(), ("C2H3".into(), "K<R".into()));
        table.insert("TEST".into(), ("N2O2-CH3".into(), "KLR".into()));
        let mut written = Vec::new();
        PTMXMLFile::new().write(&mut written, &table).unwrap();
        let back = PTMXMLFile::new().load_reader(written.as_slice()).unwrap();
        assert_eq!(back, table);
    }

    #[test]
    fn missing_fields_reuse_the_previous_ones() {
        let xml = "<PTMs><PTM><name> a </name><composition>H</composition>\
                   <possible_amino_acids>K</possible_amino_acids></PTM>\
                   <PTM><possible_amino_acids>R</possible_amino_acids></PTM></PTMs>";
        let table = PTMXMLFile::new().load_reader(xml.as_bytes()).unwrap();
        assert_eq!(table.len(), 1);
        assert_eq!(table["a"], ("H".to_owned(), "R".to_owned()));
    }

    #[test]
    fn malformed_xml_is_a_parse_error() {
        let error = PTMXMLFile::new()
            .load_reader("<PTMs><PTM></PTMs>".as_bytes())
            .unwrap_err();
        assert!(matches!(error, Error::Parse { .. }), "{error:?}");
    }
}
