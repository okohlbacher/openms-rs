// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Unimod XML modification files: `FORMAT/UnimodXMLFile.h`, its handler
//! `FORMAT/HANDLERS/UnimodXMLHandler.h`, and `CHEMISTRY/UnimodXMLDataProvider.h`.
//!
//! The source reads `unimod.xml` (and OpenMS's `custom_mods.xml`, which uses the
//! same schema) into one `ResidueModification` per allowed site of every
//! `<umod:mod>`. [`UnimodXmlFile::load`] returns those records as
//! [`ModificationRecord`]s with an empty full identifier, as the source's
//! `UnimodXMLFile::load` leaves it unset; [`UnimodXmlDataProvider`] derives the
//! full identifier and freezes them into [`ResidueModification`]s, as the
//! source provider's `setFullId()` does. A caller-owned
//! [`crate::chemistry::ModificationsDB`] takes them through
//! `ModificationsDB::extend_records`; the global registry is the compiled table
//! generated from the same files by `tools/generate_modifications.py`.
//!
//! # API mapping
//!
//! | Source (core `bc9cc12`) | Rust |
//! |---|---|
//! | `UnimodXMLFile()`, `~UnimodXMLFile()` | [`UnimodXmlFile`] (a unit struct) |
//! | `UnimodXMLFile::load(filename, std::vector<ResidueModification*>&)` | [`UnimodXmlFile::load`] (path), [`UnimodXmlFile::load_from_reader`] |
//! | `Internal::UnimodXMLHandler` (`onStartElement`, `onEndElement`, `onCharacters`) | the private event handler of this module |
//! | `UnimodXMLDataProvider(const std::string& filename)` | [`UnimodXmlDataProvider::new`] |
//! | `UnimodXMLDataProvider::loadModifications()` | [`UnimodXmlDataProvider::load_modifications`] |
//!
//! # Preserved source conventions
//!
//! - Tags are matched by their qualified name, `umod:X` or unprefixed `X`, as
//!   the handler compares them; other prefixes are ignored.
//! - One record per `<umod:specificity>`: the site's first character as the
//!   origin, `X` for a multi-character site (`N-term`, `C-term`), the position
//!   mapped to the term specificity, an unknown position warned about and read
//!   as `Anywhere`.
//! - The delta's `mono_mass` and `avge_mass` are the difference masses; its
//!   `<umod:element>` children, `number` atoms of `symbol` (digits in the symbol
//!   are the isotope, `13C` is `(13)C`), sum to the difference formula.
//!   Elements anywhere else (`<umod:Ignore>`, `<umod:aa>`, `<umod:brick>`) are
//!   skipped.
//! - A `<umod:NeutralLoss>` with elements becomes a neutral loss of its
//!   specificity, with its `mono_mass` and `avge_mass`; one without elements
//!   (the `0` loss) is dropped.
//! - Every record of a `<umod:mod>` carries the `classification` of that mod's
//!   **last** specificity: the handler sets it on the shared modification at
//!   each specificity (its own `TODO` says so), and the copies are made at the
//!   end. The compiled table keeps each specificity's own classification.
//! - The absolute masses stay zero and the provenance is `Cv`.
//!
//! # Native differences
//!
//! - `hidden` is read from each specificity (absent is `false`); the source's
//!   `ResidueModification` has no such field. The compiled table reads it too.
//! - `File::find` resolution of a bare file name against the OpenMS data path
//!   is not done: the path is opened as given.
//! - Input is limited to [`MAX_UNIMOD_XML_BYTES`] and
//!   [`MAX_UNIMOD_RECORDS`] records, checked before the records are built.
//! - A missing required attribute, an unparsable number or formula, or a
//!   malformed document is an error; the source reports a parse error too.
//!
//! See `docs/UNIMOD_XML_SUPPORT.md`.

use crate::chemistry::{
    EmpiricalFormula, ModificationProvenance, ModificationRecord, NeutralLoss, ResidueModification,
    TermSpecificity,
};
use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Largest Unimod XML document read, in bytes (the pinned `unimod.xml` is
/// about 2.5 MB).
pub const MAX_UNIMOD_XML_BYTES: usize = 64 * 1024 * 1024;

/// Largest number of site records one document may produce (the pinned
/// `unimod.xml` gives about 3,000).
pub const MAX_UNIMOD_RECORDS: usize = 1_000_000;

/// The records and warnings of one Unimod XML document.
#[derive(Clone, Debug, Default)]
pub struct UnimodXmlLoad {
    /// One record per site, in document order, with an empty `full_id`.
    pub records: Vec<ModificationRecord>,
    /// The handler's warnings, such as an unknown specificity position.
    pub warnings: Vec<String>,
}

/// Reader for Unimod XML files: source `UnimodXMLFile`.
#[derive(Clone, Copy, Debug, Default)]
pub struct UnimodXmlFile;

impl UnimodXmlFile {
    /// Reads the file at `path`: source `UnimodXMLFile::load(filename, mods)`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when the file cannot be read (the source's
    /// `FileNotFound`), and the errors of [`Self::load_from_reader`].
    pub fn load(&self, path: impl AsRef<Path>) -> Result<UnimodXmlLoad> {
        let file = File::open(path.as_ref())?;
        self.load_from_reader(BufReader::new(file))
    }

    /// Reads a Unimod XML document from `input`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Parse`] for malformed XML, a missing required attribute
    /// (`title`, `record_id`, `site`, `position`, the delta and neutral-loss
    /// masses, an element's `symbol`), an unparsable number or element, or a
    /// record id outside `u32`, and [`Error::InvalidValue`] when the document
    /// exceeds [`MAX_UNIMOD_XML_BYTES`] or [`MAX_UNIMOD_RECORDS`].
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::format::unimod_xml::UnimodXmlFile;
    ///
    /// let xml = r#"<umod:unimod xmlns:umod="http://www.unimod.org/xmlns/schema/unimod_2">
    ///   <umod:modifications>
    ///     <umod:mod title="Oxidation" full_name="Oxidation or Hydroxylation" record_id="35">
    ///       <umod:specificity hidden="0" site="M" position="Anywhere" classification="Post-translational"/>
    ///       <umod:delta mono_mass="15.994915" avge_mass="15.9994" composition="O">
    ///         <umod:element symbol="O" number="1"/>
    ///       </umod:delta>
    ///     </umod:mod>
    ///   </umod:modifications>
    /// </umod:unimod>"#;
    /// let load = UnimodXmlFile.load_from_reader(xml.as_bytes()).unwrap();
    /// assert_eq!(load.records.len(), 1);
    /// assert_eq!(load.records[0].name, "Oxidation");
    /// assert_eq!(load.records[0].origin, Some('M'));
    /// assert_eq!(load.records[0].diff_mono_mass, 15.994915);
    /// ```
    pub fn load_from_reader(&self, input: impl BufRead) -> Result<UnimodXmlLoad> {
        let bytes = read_bounded(input)?;
        let mut reader = Reader::from_reader(bytes.as_slice());
        reader.config_mut().check_end_names = true;
        let mut handler = Handler::default();
        let mut buffer = Vec::new();
        loop {
            let event = reader
                .read_event_into(&mut buffer)
                .map_err(|e| parse_error(format!("Unimod XML: {e}")))?;
            match event {
                Event::Start(element) => handler.start(&element)?,
                Event::Empty(element) => {
                    handler.start(&element)?;
                    handler.end(element.name().as_ref())?;
                }
                Event::End(element) => handler.end(element.name().as_ref())?,
                Event::Eof => break,
                _ => {}
            }
            buffer.clear();
        }
        Ok(UnimodXmlLoad {
            records: handler.records,
            warnings: handler.warnings,
        })
    }
}

/// Modification data from a Unimod XML file: source `UnimodXMLDataProvider`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnimodXmlDataProvider {
    filename: PathBuf,
}

impl UnimodXmlDataProvider {
    /// A provider reading `filename`: source
    /// `UnimodXMLDataProvider(const std::string& filename)`. Nothing is read
    /// until [`Self::load_modifications`].
    pub fn new(filename: impl Into<PathBuf>) -> Self {
        Self {
            filename: filename.into(),
        }
    }

    /// The file this provider reads.
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// Reads the file and returns one modification per site, each with its
    /// full identifier derived (`Oxidation (M)`, `Acetyl (N-term)`,
    /// `Amidated (Protein C-term)`): source `loadModifications()`.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`UnimodXmlFile::load`], and
    /// [`Error::InvalidValue`] when a record fails
    /// [`ResidueModification::from_record`] (an empty title, a control
    /// character, an origin that is not an uppercase letter).
    pub fn load_modifications(&self) -> Result<Vec<ResidueModification>> {
        let load = UnimodXmlFile.load(&self.filename)?;
        records_to_modifications(load.records)
    }
}

/// Freezes loaded records as the provider does, deriving each full identifier.
///
/// # Errors
///
/// Returns the first error of [`ResidueModification::from_record`].
pub fn records_to_modifications(
    records: Vec<ModificationRecord>,
) -> Result<Vec<ResidueModification>> {
    records
        .into_iter()
        .map(ResidueModification::from_record)
        .collect()
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn read_bounded(input: impl BufRead) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    input
        .take(MAX_UNIMOD_XML_BYTES as u64 + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() > MAX_UNIMOD_XML_BYTES {
        return Err(Error::InvalidValue(format!(
            "Unimod XML exceeds {MAX_UNIMOD_XML_BYTES} bytes"
        )));
    }
    Ok(bytes)
}

/// Whether a qualified tag name is `umod:<local>` or `<local>`, as the source
/// handler compares them.
fn is_tag(name: &[u8], local: &str) -> bool {
    name == local.as_bytes()
        || (name.len() == local.len() + 5
            && name.starts_with(b"umod:")
            && &name[5..] == local.as_bytes())
}

/// An attribute's value after the XML parser's normalisation: literal tab,
/// newline and carriage return become spaces, then entities are expanded.
fn attribute(element: &BytesStart<'_>, key: &str) -> Result<Option<String>> {
    for attr in element.attributes() {
        let attr = attr.map_err(|e| parse_error(format!("Unimod XML attribute: {e}")))?;
        if attr.key.as_ref() == key.as_bytes() {
            let raw = std::str::from_utf8(&attr.value)
                .map_err(|_| parse_error("Unimod XML attribute is not UTF-8"))?;
            let normalized: String = raw
                .chars()
                .map(|c| {
                    if matches!(c, '\t' | '\n' | '\r') {
                        ' '
                    } else {
                        c
                    }
                })
                .collect();
            let value = quick_xml::escape::unescape(&normalized)
                .map_err(|e| parse_error(format!("Unimod XML attribute: {e}")))?;
            return Ok(Some(value.into_owned()));
        }
    }
    Ok(None)
}

fn required(element: &BytesStart<'_>, key: &str) -> Result<String> {
    attribute(element, key)?.ok_or_else(|| {
        parse_error(format!(
            "Unimod XML: <{}> lacks the attribute {key}",
            String::from_utf8_lossy(element.name().as_ref())
        ))
    })
}

fn number(element: &BytesStart<'_>, key: &str) -> Result<f64> {
    let text = required(element, key)?;
    text.trim()
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| {
            parse_error(format!(
                "Unimod XML: {key}=\"{text}\" is not a finite number"
            ))
        })
}

/// The state of `UnimodXMLHandler` between events.
#[derive(Default)]
struct Handler {
    records: Vec<ModificationRecord>,
    warnings: Vec<String>,
    /// The `<umod:mod>` being read: title, full name, record id, classification.
    modification: Option<ModificationRecord>,
    sites: Vec<char>,
    hidden: Vec<bool>,
    term_specs: Vec<TermSpecificity>,
    diff_formula: EmpiricalFormula,
    diff_mono_mass: f64,
    diff_average_mass: f64,
    collect_elements: bool,
    loss_mono_mass: f64,
    loss_average_mass: f64,
    /// The neutral losses of the specificity being read.
    losses: Vec<NeutralLoss>,
    /// The neutral losses of every specificity read so far, in order.
    losses_per_specificity: Vec<Vec<NeutralLoss>>,
}

impl Handler {
    fn start(&mut self, element: &BytesStart<'_>) -> Result<()> {
        let name = element.name();
        let tag = name.as_ref();
        if is_tag(tag, "mod") {
            self.sites.clear();
            self.hidden.clear();
            let record_id = required(element, "record_id")?;
            let record_id = record_id.trim().parse::<u32>().map_err(|_| {
                parse_error(format!(
                    "Unimod XML: record_id=\"{record_id}\" is not a record id"
                ))
            })?;
            self.modification = Some(ModificationRecord {
                provenance: ModificationProvenance::Cv,
                record_id: Some(record_id),
                name: required(element, "title")?,
                full_name: attribute(element, "full_name")?.unwrap_or_default(),
                classification: String::new(),
                ..ModificationRecord::default()
            });
            return Ok(());
        }
        if is_tag(tag, "specificity") {
            self.losses.clear();
            let classification = attribute(element, "classification")?.unwrap_or_default();
            let modification = self
                .modification
                .as_mut()
                .ok_or_else(|| parse_error("Unimod XML: <specificity> outside <mod>"))?;
            modification.classification = classification;
            let site = required(element, "site")?;
            let position = required(element, "position")?;
            let term = match position.as_str() {
                "Anywhere" => TermSpecificity::Anywhere,
                "Protein N-term" => TermSpecificity::ProteinNTerm,
                "Protein C-term" => TermSpecificity::ProteinCTerm,
                "Any C-term" => TermSpecificity::CTerm,
                "Any N-term" => TermSpecificity::NTerm,
                other => {
                    self.warnings.push(format!(
                        "Don't know allowed position called: '{other}' - setting to anywhere"
                    ));
                    TermSpecificity::Anywhere
                }
            };
            self.term_specs.push(term);
            let mut chars = site.chars();
            let origin = match (chars.next(), chars.next()) {
                (Some(only), None) => only,
                (Some(_), Some(_)) => 'X',
                (None, _) => return Err(parse_error("Unimod XML: empty specificity site")),
            };
            self.sites.push(origin);
            self.hidden
                .push(attribute(element, "hidden")?.is_some_and(|flag| flag.trim() == "1"));
            if self.sites.len() > MAX_UNIMOD_RECORDS
                || self.records.len().saturating_add(self.sites.len()) > MAX_UNIMOD_RECORDS
            {
                return Err(Error::InvalidValue(format!(
                    "Unimod XML has more than {MAX_UNIMOD_RECORDS} site records"
                )));
            }
            return Ok(());
        }
        if is_tag(tag, "NeutralLoss") {
            self.loss_mono_mass = number(element, "mono_mass")?;
            self.loss_average_mass = number(element, "avge_mass")?;
            self.collect_elements = true;
            return Ok(());
        }
        if is_tag(tag, "delta") {
            self.diff_average_mass = number(element, "avge_mass")?;
            self.diff_mono_mass = number(element, "mono_mass")?;
            self.collect_elements = true;
            return Ok(());
        }
        if is_tag(tag, "element") && self.collect_elements {
            let symbol = required(element, "symbol")?;
            let count = required(element, "number")?;
            let isotope: String = symbol.chars().filter(char::is_ascii_digit).collect();
            let bare: String = symbol.chars().filter(|c| !c.is_ascii_digit()).collect();
            let text = if isotope.is_empty() {
                format!("{bare}{count}")
            } else {
                format!("({isotope}){bare}{count}")
            };
            let atoms = EmpiricalFormula::parse(&text)
                .map_err(|e| parse_error(format!("Unimod XML: element {symbol} x {count}: {e}")))?;
            self.diff_formula = self.diff_formula.checked_add(&atoms)?;
        }
        Ok(())
    }

    fn end(&mut self, tag: &[u8]) -> Result<()> {
        if is_tag(tag, "mod") {
            let modification = self
                .modification
                .take()
                .ok_or_else(|| parse_error("Unimod XML: unbalanced </mod>"))?;
            let formula = std::mem::take(&mut self.diff_formula);
            let losses = std::mem::take(&mut self.losses_per_specificity);
            for (index, origin) in self.sites.iter().enumerate() {
                let mut record = modification.clone();
                record.diff_mono_mass = self.diff_mono_mass;
                record.diff_average_mass = self.diff_average_mass;
                record.diff_formula = formula.clone();
                record.origin = Some(*origin);
                record.term_specificity = self.term_specs[index];
                record.hidden = self.hidden[index];
                record.neutral_losses = losses.get(index).cloned().unwrap_or_default();
                self.records.push(record);
            }
            self.diff_mono_mass = 0.0;
            self.diff_average_mass = 0.0;
            self.collect_elements = false;
            self.term_specs.clear();
            self.sites.clear();
            self.hidden.clear();
            return Ok(());
        }
        if is_tag(tag, "specificity") {
            self.losses_per_specificity
                .push(std::mem::take(&mut self.losses));
            return Ok(());
        }
        if is_tag(tag, "delta") {
            self.collect_elements = false;
            return Ok(());
        }
        if is_tag(tag, "NeutralLoss") {
            self.collect_elements = false;
            if !self.diff_formula.is_empty() {
                let formula = std::mem::take(&mut self.diff_formula);
                self.losses.push(NeutralLoss::new(
                    formula,
                    self.loss_mono_mass,
                    self.loss_average_mass,
                )?);
            }
        }
        Ok(())
    }
}
