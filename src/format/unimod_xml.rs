// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Unimod XML: the modification database of unimod.org, read at run time.
//!
//! Ports `FORMAT/UnimodXMLFile.h`, its handler `FORMAT/HANDLERS/UnimodXMLHandler.h`
//! and `CHEMISTRY/UnimodXMLDataProvider.h`, the provider that wraps the file
//! for a modification registry.
//!
//! The default build does not need this: `ModificationsDB::global()` is built
//! from a table `tools/generate_modifications.py` flattens from the pinned
//! `resources/modifications/unimod.xml` offline. This reader produces the same
//! records from any Unimod-style file at run time, so a caller can load a newer
//! `unimod.xml`, the OpenMS `custom_mods.xml`, or its own file, and build a
//! registry with [`ModificationsDB::from_records`].
//!
//! One `<umod:mod>` becomes one [`ResidueModification`] per
//! `<umod:specificity>`, in document order. See
//! `docs/UNIMOD_XML_SUPPORT.md` for the member mapping and evidence.
//!
//! [`ModificationsDB::from_records`]: crate::chemistry::ModificationsDB::from_records
//! [`ResidueModification`]: crate::chemistry::ResidueModification

use super::sax_xml::{SaxAttributes, SaxEvent, SaxLimits, parse};
use crate::chemistry::{
    EmpiricalFormula, ModificationProvenance, ModificationRecord, NeutralLoss, ResidueModification,
    TermSpecificity,
};
use crate::{Error, Result};
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Resource bounds for reading one Unimod XML document. Exceeding any of them
/// is [`Error::InvalidValue`] and yields no result.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Decompressed input bytes. The pinned `unimod.xml` is 2.4 MB.
    pub max_input_bytes: u64,
    /// Open elements at any one time.
    pub max_depth: usize,
    /// Parser events in the document.
    pub max_events: u64,
    /// Modification records in the result (one per specificity).
    pub max_records: usize,
    /// Specificities of one `<umod:mod>`.
    pub max_specificities: usize,
    /// Neutral losses of one specificity.
    pub max_neutral_losses: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_input_bytes: 256 * 1024 * 1024,
            max_depth: 64,
            max_events: 50_000_000,
            max_records: 1_000_000,
            max_specificities: 10_000,
            max_neutral_losses: 10_000,
        }
    }
}

/// What to do with the two specificity attributes the source handler drops.
///
/// Unimod gives every `<umod:specificity>` its own `classification` and
/// `hidden` flag. Source `UnimodXMLHandler` keeps neither per specificity: it
/// ignores `hidden`, and it writes each specificity's classification into the
/// one modification all of them are copied from, so the *last* one wins for
/// every site (`UnimodXMLHandler.cpp:60-63`, marked `TODO` there).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SpecificityAttributes {
    /// Keep each specificity's own `classification` text and `hidden` flag.
    ///
    /// This loses nothing, and it is what the embedded table behind
    /// `ModificationsDB::global()` records, so reading the pinned
    /// `unimod.xml` this way reproduces that registry's records exactly.
    #[default]
    Retain,
    /// The source's result: `hidden` is `false`, and every specificity gets
    /// the classification of the modification's last specificity, as source
    /// `ResidueModification::setSourceClassification` canonicalises it —
    /// matched case-insensitively against the Unimod names and returned in
    /// `getSourceClassificationName`'s spelling, or the empty string for a
    /// name it does not know.
    Source,
}

/// Options for reading Unimod XML.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReadOptions {
    /// How specificity classifications and hidden flags are carried.
    pub specificity_attributes: SpecificityAttributes,
    /// Resource bounds.
    pub limits: ReadLimits,
}

impl ReadOptions {
    /// Options that reproduce the source handler's records.
    pub fn source() -> Self {
        Self {
            specificity_attributes: SpecificityAttributes::Source,
            limits: ReadLimits::default(),
        }
    }

    fn sax(&self) -> SaxLimits {
        SaxLimits {
            max_input_bytes: self.limits.max_input_bytes,
            max_depth: self.limits.max_depth,
            max_events: self.limits.max_events,
            max_attributes: 256,
            max_text_bytes: 16 * 1024 * 1024,
        }
    }
}

/// Source `ResidueModification::setSourceClassification(const String&)`
/// followed by `getSourceClassificationName()`.
fn source_classification(text: &str) -> &'static str {
    match text.to_ascii_lowercase().as_str() {
        "artifact" | "artefact" => "Artefact",
        "natural" => "Natural",
        "hypothetical" => "Hypothetical",
        "post-translational" => "Post-translational",
        "multiple" => "Multiple",
        "chemical derivative" => "Chemical derivative",
        "isotopic label" => "Isotopic label",
        "pre-translational" => "Pre-translational",
        "other glycosylation" => "Other glycosylation",
        "n-linked glycosylation" => "N-linked glycosylation",
        "aa substitution" => "AA substitution",
        "other" => "Other",
        "non-standard residue" => "Non-standard residue",
        "co-translational" => "Co-translational",
        "o-linked glycosylation" => "O-linked glycosylation",
        _ => "",
    }
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// Source `attributeAsString_`: a missing attribute is a fatal parse error.
fn required<'a>(attributes: &'a SaxAttributes, name: &str) -> Result<&'a str> {
    attributes.required(name)
}

/// A finite decimal through the source's `StringUtils::toDouble`
/// (`attributeAsDouble_`, and the direct calls of `UnimodXMLHandler.cpp:122-123`);
/// it also accepts `nan` and `inf`, which a modification mass cannot be, so
/// those are refused.
fn mass(attributes: &SaxAttributes, name: &str) -> Result<f64> {
    let text = required(attributes, name)?;
    crate::data_structures::string_utils::to_double(text)
        .ok()
        .filter(|value| value.is_finite())
        .ok_or_else(|| {
            parse_error(format!(
                "attribute '{name}' is not a finite number: {text:?}"
            ))
        })
}

fn boolean(text: &str) -> Result<bool> {
    match text.trim_matches([' ', '\t', '\n', '\r']) {
        "1" | "true" => Ok(true),
        "0" | "false" => Ok(false),
        _ => Err(parse_error(format!(
            "attribute 'hidden' is not a boolean: {text:?}"
        ))),
    }
}

/// Everything the source keeps per `<umod:specificity>`.
struct Specificity {
    origin: char,
    term: TermSpecificity,
    classification: String,
    hidden: bool,
    neutral_losses: Vec<NeutralLoss>,
}

/// The modification-level fields the source stores in its `modification_`.
#[derive(Default)]
struct Modification {
    title: String,
    full_name: String,
    record_id: u32,
    /// The classification of the specificity seen last, as the source's
    /// single `setSourceClassification` call leaves it.
    last_classification: String,
}

/// The state of source `UnimodXMLHandler`, member for member where the
/// member affects the result.
#[derive(Default)]
struct Handler {
    options: ReadOptions,
    /// `modification_`: `None` outside a `<umod:mod>`.
    modification: Option<Modification>,
    /// `avge_mass_` / `mono_mass_`: the last `<umod:delta>` masses.
    average_mass: f64,
    mono_mass: f64,
    /// `diff_formula_`: element counts collected while `collect_elements`.
    diff_formula: EmpiricalFormula,
    /// `collect_elements_`: inside `<umod:delta>` or `<umod:NeutralLoss>`.
    collect_elements: bool,
    /// Masses of the open `<umod:NeutralLoss>`.
    loss_mono_mass: f64,
    loss_average_mass: f64,
    /// Neutral losses of the open specificity (`neutral_loss_diff_formula_`
    /// with the two parallel mass vectors).
    losses: Vec<NeutralLoss>,
    /// Closed specificities of the open modification (`sites_`,
    /// `term_specs_` and the per-specificity loss vectors). The open one is
    /// `open_specificity`, which the source pushes at its start tag.
    specificities: Vec<Specificity>,
    open_specificity: Option<Specificity>,
    output: Vec<ResidueModification>,
}

/// The local part of a tag the source matches in both spellings, with and
/// without the `umod:` prefix, and only those two.
fn unimod_tag(qname: &str) -> &str {
    qname.strip_prefix("umod:").unwrap_or(qname)
}

impl Handler {
    fn event(&mut self, event: SaxEvent<'_>) -> Result<()> {
        match event {
            SaxEvent::Start(qname, attributes) => self.start(qname, attributes),
            SaxEvent::End(qname) => self.end(qname),
            SaxEvent::Characters(_) => Ok(()),
        }
    }

    fn start(&mut self, qname: &str, attributes: &SaxAttributes) -> Result<()> {
        match unimod_tag(qname) {
            "mod" => {
                // The source discards (without freeing) a modification still
                // open here; a nested <mod> cannot close the outer one twice.
                if self.modification.is_some() {
                    return Err(parse_error("<umod:mod> nested in <umod:mod>"));
                }
                self.specificities.clear();
                let record_id = required(attributes, "record_id")?;
                let record_id = record_id
                    .trim_matches([' ', '\t', '\n', '\r'])
                    .parse::<u32>()
                    .map_err(|_| {
                        parse_error(format!(
                            "attribute 'record_id' is not a record number: {record_id:?}"
                        ))
                    })?;
                self.modification = Some(Modification {
                    title: required(attributes, "title")?.to_owned(),
                    full_name: required(attributes, "full_name")?.to_owned(),
                    record_id,
                    last_classification: String::new(),
                });
            }
            "specificity" => {
                let modification = self
                    .modification
                    .as_mut()
                    .ok_or_else(|| parse_error("<umod:specificity> outside <umod:mod>"))?;
                if self.open_specificity.is_some() {
                    return Err(parse_error(
                        "<umod:specificity> nested in <umod:specificity>",
                    ));
                }
                if self.specificities.len() == self.options.limits.max_specificities {
                    return Err(Error::InvalidValue(
                        "Unimod specificity limit exceeded".into(),
                    ));
                }
                self.losses.clear();
                let classification = required(attributes, "classification")?;
                modification.last_classification = classification.to_owned();
                let site = required(attributes, "site")?;
                let position = required(attributes, "position")?;
                let term = match position {
                    "Anywhere" => TermSpecificity::Anywhere,
                    "Protein N-term" => TermSpecificity::ProteinNTerm,
                    "Protein C-term" => TermSpecificity::ProteinCTerm,
                    "Any C-term" => TermSpecificity::CTerm,
                    "Any N-term" => TermSpecificity::NTerm,
                    // Source: a warning "Don't know allowed position called:
                    // '...' - setting to anywhere", and the record is kept.
                    _ => TermSpecificity::Anywhere,
                };
                // Source: a site longer than one character ("N-term",
                // "C-term") is the wildcard X; otherwise its first character.
                let mut characters = site.chars();
                let origin = match (characters.next(), characters.next()) {
                    (Some(c), None) => c,
                    (Some(_), Some(_)) => 'X',
                    // The source reads site[0] of an empty string: the
                    // terminating NUL, which no residue has.
                    (None, _) => return Err(parse_error("empty specificity site")),
                };
                let hidden = match self.options.specificity_attributes {
                    SpecificityAttributes::Retain => boolean(required(attributes, "hidden")?)?,
                    SpecificityAttributes::Source => false,
                };
                self.open_specificity = Some(Specificity {
                    origin,
                    term,
                    classification: classification.to_owned(),
                    hidden,
                    neutral_losses: Vec::new(),
                });
            }
            "NeutralLoss" => {
                self.loss_mono_mass = mass(attributes, "mono_mass")?;
                self.loss_average_mass = mass(attributes, "avge_mass")?;
                self.collect_elements = true;
            }
            "delta" => {
                // Source reads these through valueByIndex(index(...)) with no
                // presence check; an absent one is undefined there.
                self.average_mass = mass(attributes, "avge_mass")?;
                self.mono_mass = mass(attributes, "mono_mass")?;
                self.collect_elements = true;
            }
            "element" if self.collect_elements => {
                let symbol = required(attributes, "symbol")?;
                let number = required(attributes, "number")?;
                // Source: digits of the symbol are the isotope, the rest the
                // element, giving "(13)C6" for symbol "13C", number "6".
                let (isotope, element): (String, String) =
                    symbol.chars().partition(char::is_ascii_digit);
                let formula = if isotope.is_empty() {
                    format!("{element}{number}")
                } else {
                    format!("({isotope}){element}{number}")
                };
                let formula = EmpiricalFormula::parse(&formula).map_err(|error| {
                    parse_error(format!(
                        "invalid Unimod element {symbol:?} x {number:?}: {error}"
                    ))
                })?;
                self.diff_formula = self.diff_formula.checked_add(&formula)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn end(&mut self, qname: &str) -> Result<()> {
        match unimod_tag(qname) {
            "mod" => self.finish_modification()?,
            "specificity" => {
                // Source: the per-specificity loss vectors are pushed here.
                if let Some(mut specificity) = self.open_specificity.take() {
                    specificity.neutral_losses = std::mem::take(&mut self.losses);
                    self.specificities.push(specificity);
                }
            }
            "delta" => self.collect_elements = false,
            "NeutralLoss" => {
                self.collect_elements = false;
                // Source: a loss without element children adds nothing, and
                // the collected formula is reset only when it is kept.
                if !self.diff_formula.is_empty() {
                    if self.losses.len() == self.options.limits.max_neutral_losses {
                        return Err(Error::InvalidValue(
                            "Unimod neutral loss limit exceeded".into(),
                        ));
                    }
                    let formula = std::mem::take(&mut self.diff_formula);
                    self.losses.push(NeutralLoss::new(
                        formula,
                        self.loss_mono_mass,
                        self.loss_average_mass,
                    )?);
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn finish_modification(&mut self) -> Result<()> {
        let Some(modification) = self.modification.take() else {
            return Err(parse_error("</umod:mod> without <umod:mod>"));
        };
        let diff_formula = std::mem::take(&mut self.diff_formula);
        let count = self
            .output
            .len()
            .checked_add(self.specificities.len())
            .filter(|&n| n <= self.options.limits.max_records)
            .ok_or_else(|| Error::InvalidValue("Unimod record limit exceeded".into()))?;
        self.output.reserve(count - self.output.len());
        for specificity in self.specificities.drain(..) {
            let (classification, hidden) = match self.options.specificity_attributes {
                SpecificityAttributes::Retain => (specificity.classification, specificity.hidden),
                SpecificityAttributes::Source => (
                    source_classification(&modification.last_classification).to_owned(),
                    false,
                ),
            };
            let record = ModificationRecord {
                provenance: ModificationProvenance::Cv,
                record_id: Some(modification.record_id),
                name: modification.title.clone(),
                full_name: modification.full_name.clone(),
                full_id: String::new(),
                origin: Some(specificity.origin),
                term_specificity: specificity.term,
                diff_formula: diff_formula.clone(),
                diff_mono_mass: self.mono_mass,
                diff_average_mass: self.average_mass,
                hidden,
                classification,
                neutral_losses: specificity.neutral_losses,
                ..ModificationRecord::default()
            };
            if record.name.is_empty() {
                // Source setFullId: "Cannot create full ID for modification
                // with missing (short) ID." (ResidueModification.cpp:140).
                return Err(Error::MissingInformation(
                    "Cannot create full ID for modification with missing (short) ID.".into(),
                ));
            }
            self.output
                .push(ResidueModification::from_record(record).map_err(|error| {
                    parse_error(format!(
                        "invalid Unimod modification {:?}: {error}",
                        modification.title
                    ))
                })?);
        }
        self.average_mass = 0.0;
        self.mono_mass = 0.0;
        self.collect_elements = false;
        Ok(())
    }
}

/// Read a Unimod XML document with the default [`ReadOptions`].
///
/// Source `UnimodXMLFile::load`, which appends to a caller's vector of raw
/// pointers; this returns the records instead, and only on success. Records
/// come in document order: modifications as they appear, and within one, its
/// specificities as they appear.
///
/// Each record carries the modification's `title` as its name, `full_name`,
/// `record_id`, the `<umod:delta>` masses and element composition (as the
/// difference formula), and its specificity's site, position and neutral
/// losses. It has [`ModificationProvenance::Cv`], as the source's
/// `setProvenance(CV)`.
///
/// The reader is the source handler's state machine:
///
/// - Tags are matched with and without the `umod:` prefix, as the source's
///   namespace-unaware parse sees them; any other prefix is not recognised.
/// - Only `<umod:element>` children of `<umod:delta>` and `<umod:NeutralLoss>`
///   contribute to formulas; those of `<umod:Ignore>`, `<umod:aa>` and
///   `<umod:brick>` do not (OpenMS issue 10030).
/// - A neutral loss without element children is dropped, as the source drops
///   it; its masses come from its own `mono_mass` and `avge_mass`.
/// - A site of more than one character (Unimod's `N-term`, `C-term`) is the
///   wildcard `X`, which a terminal specificity then stores as no origin.
/// - An unknown `position` is `Anywhere`; the source logs a warning.
/// - A `<umod:delta>` outside every `<umod:mod>` still sets the masses the
///   next modification uses unless that one has its own, as in the source.
///
/// Native differences: every record has its full identifier (for example
/// `Oxidation (M)`) derived at construction, which the source leaves to
/// [`UnimodXMLDataProvider`]'s `setFullId`; masses must be finite; an XML
/// structure the source mishandles through a dangling pointer (a
/// specificity outside a modification, nested modifications) is refused.
///
/// # Errors
///
/// [`Error::Parse`] for XML that is not well formed, a missing required
/// attribute, a mass or record number that does not parse, or a formula
/// element the element table does not hold — each an `Exception::ParseError`
/// in the source; [`Error::MissingInformation`] for a modification with an
/// empty `title`; [`Error::Unsupported`] for a DTD or a non-UTF-8 encoding;
/// [`Error::InvalidValue`] when a bound in [`ReadLimits`] is exceeded;
/// [`Error::Io`] when the input fails.
pub fn read(input: impl BufRead) -> Result<Vec<ResidueModification>> {
    read_with_options(input, &ReadOptions::default())
}

/// [`read`] with explicit options.
///
/// # Errors
///
/// As [`read`].
pub fn read_with_options(
    input: impl BufRead,
    options: &ReadOptions,
) -> Result<Vec<ResidueModification>> {
    let mut handler = Handler {
        options: *options,
        ..Handler::default()
    };
    parse(input, &options.sax(), |event| handler.event(event))?;
    Ok(handler.output)
}

/// Load a Unimod XML file with the default [`ReadOptions`]; gzip and bzip2
/// input are recognised by their leading bytes.
///
/// Source `UnimodXMLFile::load` first resolves `filename` with `File::find`,
/// so `"CHEMISTRY/unimod.xml"` finds the installed copy. This reads the path
/// it is given; resolve a shared-data name with
/// `system::file::FileContext::find` first.
///
/// # Errors
///
/// [`Error::Io`] when the file cannot be opened (the source's
/// `Exception::FileNotFound`); otherwise as [`read`].
pub fn load(path: impl AsRef<Path>) -> Result<Vec<ResidueModification>> {
    load_with_options(path, &ReadOptions::default())
}

/// [`load`] with explicit options.
///
/// # Errors
///
/// As [`load`].
pub fn load_with_options(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<Vec<ResidueModification>> {
    read_with_options(super::path_io::open(path.as_ref())?, options)
}

/// A source of modification records backed by one Unimod XML file.
///
/// Source `UnimodXMLDataProvider`, a `ModificationDataProvider` that wraps
/// `UnimodXMLFile::load` and then gives every record its full identifier.
/// Construction stores the path only; nothing is read until
/// [`load_modifications`](Self::load_modifications). Feed the records to
/// [`ModificationsDB::from_records`] or
/// [`ModificationsDB::extend_records`].
///
/// [`ModificationsDB::from_records`]: crate::chemistry::ModificationsDB::from_records
/// [`ModificationsDB::extend_records`]: crate::chemistry::ModificationsDB::extend_records
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnimodXMLDataProvider {
    filename: PathBuf,
    options: ReadOptions,
}

impl UnimodXMLDataProvider {
    /// A provider for the Unimod XML file at `filename`, read with the default
    /// [`ReadOptions`].
    ///
    /// The source resolves `filename` through `File::find` when it loads;
    /// this reads the path as given (see [`load`]).
    pub fn new(filename: impl Into<PathBuf>) -> Self {
        Self::with_options(filename, ReadOptions::default())
    }

    /// A provider that reads with explicit options, for example
    /// [`ReadOptions::source`] to reproduce the source's records.
    pub fn with_options(filename: impl Into<PathBuf>, options: ReadOptions) -> Self {
        Self {
            filename: filename.into(),
            options,
        }
    }

    /// The file this provider reads.
    pub fn filename(&self) -> &Path {
        &self.filename
    }

    /// The options this provider reads with.
    pub fn options(&self) -> &ReadOptions {
        &self.options
    }

    /// Read the file and return its records, each with its full identifier.
    ///
    /// Source `loadModifications`, which returns owning pointers; the records
    /// are owned values here, and a failure returns none of them.
    ///
    /// # Errors
    ///
    /// As [`load`].
    pub fn load_modifications(&self) -> Result<Vec<ResidueModification>> {
        load_with_options(&self.filename, &self.options)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ONE: &str = r#"<umod:unimod xmlns:umod="u"><umod:modifications>
<umod:mod title="T" full_name="Test" record_id="7">
  <umod:specificity hidden="1" site="S" position="Anywhere" classification="Artefact">
    <umod:NeutralLoss mono_mass="0" avge_mass="0" flag="false" composition="0"/>
    <umod:NeutralLoss mono_mass="97.976896" avge_mass="97.9952" flag="false" composition="H(3) O(4) P">
      <umod:element symbol="H" number="3"/><umod:element symbol="O" number="4"/><umod:element symbol="P" number="1"/>
    </umod:NeutralLoss>
  </umod:specificity>
  <umod:specificity hidden="0" site="N-term" position="Protein N-term" classification="Post-translational"/>
  <umod:delta mono_mass="79.966331" avge_mass="79.9799" composition="H O(3) P">
    <umod:element symbol="H" number="1"/><umod:element symbol="O" number="3"/><umod:element symbol="P" number="1"/>
  </umod:delta>
  <umod:Ignore><umod:element symbol="C" number="100"/></umod:Ignore>
</umod:mod></umod:modifications></umod:unimod>"#;

    #[test]
    fn one_record_per_specificity_with_source_state_rules() {
        let records = read(ONE.as_bytes()).unwrap();
        assert_eq!(records.len(), 2);
        let first = &records[0];
        assert_eq!(first.name(), "T");
        assert_eq!(first.full_id(), "T (S)");
        assert_eq!(first.record_id(), Some(7));
        assert_eq!(first.diff_formula().to_string(), "H1O3P1");
        assert_eq!(first.diff_mono_mass(), 79.966331);
        assert!(first.is_hidden());
        assert_eq!(first.classification(), "Artefact");
        assert_eq!(first.neutral_losses().len(), 1);
        assert_eq!(first.neutral_losses()[0].mono_mass(), 97.976896);
        let second = &records[1];
        assert_eq!(second.origin(), None);
        assert_eq!(second.term_specificity(), TermSpecificity::ProteinNTerm);
        assert_eq!(second.full_id(), "T (Protein N-term)");
        assert!(second.neutral_losses().is_empty());
        assert_eq!(second.classification(), "Post-translational");

        let source = read_with_options(ONE.as_bytes(), &ReadOptions::source()).unwrap();
        for record in &source {
            assert!(!record.is_hidden());
            assert_eq!(record.classification(), "Post-translational");
        }
    }

    #[test]
    fn classification_names_follow_the_source_table() {
        assert_eq!(source_classification("ARTIFACT"), "Artefact");
        assert_eq!(source_classification("aa substitution"), "AA substitution");
        assert_eq!(source_classification("Something else"), "");
    }

    #[test]
    fn malformed_records_are_refused() {
        for text in [
            r#"<m><mod title="T" full_name="F"/></m>"#,
            r#"<m><specificity site="S" position="Anywhere" classification="x" hidden="0"/></m>"#,
            r#"<m><mod title="" full_name="F" record_id="1"><specificity site="S" position="Anywhere" classification="x" hidden="0"/></mod></m>"#,
            r#"<m><mod title="T" full_name="F" record_id="1"><delta mono_mass="nan" avge_mass="1"/></mod></m>"#,
        ] {
            assert!(read(text.as_bytes()).is_err(), "{text}");
        }
    }
}
