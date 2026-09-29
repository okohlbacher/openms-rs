// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! The TraML reader and writer behind [`traml`](crate::format::traml).
//!
//! Ports `FORMAT/HANDLERS/TraMLHandler.h` (with its `.cpp`). The source is a
//! Xerces SAX handler that the header marks "Do not use this class. It is
//! only needed in TraMLFile"; here it is two functions,
//! [`read_document`] and [`write_document`], that `traml` calls and that are
//! public so a caller can read from or write to any stream.
//!
//! The reader keeps the source's structure: one *pending* object per element
//! type (`actual_protein_`, `actual_transition_`, ...), filled while its
//! element is open and appended to the experiment when it closes, with every
//! `cvParam` and `userParam` routed by its parent and grandparent tag exactly
//! as `handleCVParam_` and `handleUserParam_` route it.
//!
//! Where the source logs a problem and drops data (an unknown element, a
//! `cvParam` or `userParam` in a place it does not handle, a `RetentionTime`,
//! `Configuration` or `Target` under an unexpected parent), this port refuses
//! the document by default. [`ReadOptions::source_compatible`] selects the
//! source behaviour instead: the data is dropped and the message the source
//! logs is returned in [`TraMLLoad::warnings`]. Messages the source logs
//! without dropping anything (an obsolete CV term, a CV term name that does
//! not match the vocabulary, an unsupported retention-time term) are always
//! returned as warnings.
//!
//! See [the support document](https://github.com/OpenMS/OpenMS4-R/blob/main/docs/TRAML_SUPPORT.md)
//! for the member mapping and the evidence. The source contains no
//! `#pragma omp`; the port is serial.

use crate::analysis::openswath::include_exclude_target::IncludeExcludeTarget;
use crate::analysis::openswath::reaction_monitoring_transition::{
    DecoyTransitionType, ReactionMonitoringTransition,
};
use crate::analysis::openswath::targeted_experiment::TargetedExperiment;
use crate::analysis::openswath::targeted_experiment_helper::{
    CV, Compound, Configuration, Contact, Instrument, Interpretation, Peptide, PeptideModification,
    Prediction, Protein, Publication, RTType, RTUnit, RetentionTime, TraMLProduct,
};
use crate::chemistry::{ModificationsDB, PeptideFragmentType, TermSpecificity};
use crate::format::controlled_vocabulary::{ControlledVocabulary, XRefType};
use crate::format::file_info::text_format::ostream_g;
use crate::metadata::{
    CVTerm, CVTermList, MetaInfo, MetaValue, MetaValueData, Software, SourceFile, Unit,
};
use crate::{Error, Result};
use quick_xml::Reader;
use quick_xml::events::{BytesStart, Event};
use std::io::{BufRead, Read, Write};

/// Ceilings and the compatibility switch for [`read_document`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadOptions {
    /// Largest accepted (decompressed) document, in bytes. Default 4 GiB.
    pub max_input_bytes: u64,
    /// Largest number of XML elements. Default 100 million.
    pub max_elements: usize,
    /// Deepest element nesting. Default 64; TraML nests at most 8 deep.
    pub max_depth: usize,
    /// Reproduce the source's log-and-continue behaviour: data the source
    /// drops is dropped and reported in [`TraMLLoad::warnings`] instead of
    /// refusing the document. Default `false`.
    pub source_compatible: bool,
}

impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            max_input_bytes: 4 * 1024 * 1024 * 1024,
            max_elements: 100_000_000,
            max_depth: 64,
            source_compatible: false,
        }
    }
}

impl ReadOptions {
    /// The defaults with [`source_compatible`](Self::source_compatible) set,
    /// which is what the source `TraMLFile::load` does.
    pub fn source() -> Self {
        Self {
            source_compatible: true,
            ..Self::default()
        }
    }
}

/// A loaded experiment and the messages the source would have logged.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TraMLLoad {
    /// The experiment read from the document.
    pub experiment: TargetedExperiment,
    /// The source's `warning(LOAD, ...)` and `error(LOAD, ...)` messages, in
    /// document order, without the `While loading '<file>': ` prefix.
    pub warnings: Vec<String>,
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// `StringUtils::toDouble`: surrounding whitespace is ignored; `nan` and
/// `inf` are accepted, as `std::from_chars` accepts them.
pub(crate) fn source_double(text: &str) -> Result<f64> {
    let trimmed = text.trim_matches([' ', '\t', '\n', '\r']);
    trimmed
        .parse::<f64>()
        .map_err(|_| parse_error(format!("Could not convert string '{text}' to a double value")))
}

/// `StringUtils::toInt32`: surrounding whitespace is ignored and the whole
/// remainder must be a decimal integer in the `i32` range.
pub(crate) fn source_int32(text: &str) -> Result<i32> {
    let trimmed = text.trim_matches([' ', '\t', '\n', '\r']);
    trimmed
        .parse::<i32>()
        .map_err(|_| parse_error(format!("Could not convert string '{text}' to an int value")))
}

/// `StringUtils::toInt64`, as [`source_int32`] for `i64`.
fn source_int64(text: &str) -> Result<i64> {
    let trimmed = text.trim_matches([' ', '\t', '\n', '\r']);
    trimmed
        .parse::<i64>()
        .map_err(|_| parse_error(format!("Could not convert string '{text}' to an int64 value")))
}

/// `XMLHandler::fromXSDString`: the typed value of a `userParam`.
///
/// `xsd:double`, `xsd:float` and `xsd:decimal` become a float; `xsd:byte`,
/// `xsd:int`, `xsd:short` and their unsigned forms an integer read as 32
/// bits; `xsd:long`, `xsd:unsignedLong`, `xsd:integer` and the four signed
/// `*Integer` types an integer read as 64 bits; anything else a string. A
/// non-finite float is stored as the source stores it.
///
/// # Errors
///
/// [`Error::Parse`] when a numeric type's value does not convert, where the
/// source throws `Exception::ConversionError`.
pub fn from_xsd_string(xsd_type: &str, value: &str) -> Result<MetaValue> {
    Ok(match xsd_type {
        "xsd:double" | "xsd:float" | "xsd:decimal" => MetaValue::source_float(source_double(value)?),
        "xsd:byte" | "xsd:int" | "xsd:unsignedShort" | "xsd:short" | "xsd:unsignedByte"
        | "xsd:unsignedInt" => MetaValue::from(source_int32(value)?),
        "xsd:long"
        | "xsd:unsignedLong"
        | "xsd:integer"
        | "xsd:negativeInteger"
        | "xsd:nonNegativeInteger"
        | "xsd:nonPositiveInteger"
        | "xsd:positiveInteger" => MetaValue::from(source_int64(value)?),
        _ => MetaValue::from(value),
    })
}

/// `DataValue::toString()` (full precision): a string verbatim, an integer
/// in decimal, a float as `StringUtils::toStr(double)`, a list as `[a, b]`
/// with each element converted likewise, and an empty value as the empty
/// string. A unit is not part of the text, as in the source.
pub(crate) fn data_value_text(value: &MetaValue) -> String {
    use crate::param::value::format_float;
    fn list<T>(values: &[T], element: impl Fn(&T) -> String) -> String {
        let parts: Vec<String> = values.iter().map(element).collect();
        format!("[{}]", parts.join(", "))
    }
    match value.data() {
        MetaValueData::Empty => String::new(),
        MetaValueData::String(text) => text.clone(),
        MetaValueData::Integer(number) => number.to_string(),
        MetaValueData::Float(number) => format_float(*number, true),
        MetaValueData::StringList(values) => list(values, Clone::clone),
        MetaValueData::IntegerList(values) => list(values, ToString::to_string),
        MetaValueData::FloatList(values) => list(values, |x| format_float(*x, true)),
    }
}

/// `XMLHandler::writeXMLEscape`: `&`, `>`, `"`, `<` and `'` as entities.
pub fn write_xml_escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '<' => out.push_str("&lt;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

/// Escaping for the places the source writes a string without
/// `writeXMLEscape`: only the three characters that would make the document
/// malformed or change its meaning (`&`, `<`, `"`), so every document the
/// source writes well-formed is written byte-identical.
fn attribute_minimal_escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '"' => out.push_str("&quot;"),
            _ => out.push(c),
        }
    }
    out
}

/// The five source tags `onStartElement` skips before anything else.
const START_IGNORED: &[&str] = &[
    "TraML",
    "ContactList",
    "CompoundList",
    "TransitionList",
    "ConfigurationList",
    "cvList",
    "InstrumentList",
    "SoftwareList",
    "PublicationList",
    "ProteinList",
    "SourceFileList",
    "InterpretationList",
    "Evidence",
    "ValidationStatus",
    "Sequence",
    "Precursor",
    "Product",
    "IntermediateProduct",
    "TargetIncludeList",
    "TargetExcludeList",
    "TargetList",
    "RetentionTimeList",
];

/// The tags `onEndElement` skips.
const END_IGNORED: &[&str] = &[
    "TraML",
    "ContactList",
    "CompoundList",
    "TransitionList",
    "ConfigurationList",
    "cvList",
    "InstrumentList",
    "SoftwareList",
    "PublicationList",
    "ProteinList",
    "SourceFileList",
    "InterpretationList",
    "Evidence",
    "cvParam",
    "userParam",
    "cv",
    "Sequence",
    "Precursor",
    "RetentionTimeList",
    "TargetList",
    "TargetIncludeList",
    "TargetExcludeList",
    "ProteinRef",
    "Modification",
];

/// The attributes of one start tag, unescaped and whitespace-normalised.
struct Attributes(Vec<(String, String)>);

impl Attributes {
    fn optional(&self, name: &str) -> Option<&str> {
        self.0
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.as_str())
    }

    /// `attributeAsString_`: a missing attribute is fatal.
    fn required(&self, name: &str) -> Result<&str> {
        self.optional(name)
            .ok_or_else(|| parse_error(format!("Required attribute '{name}' not present!")))
    }
}

/// The reader's state: the source handler's `actual_*_` members.
#[derive(Default)]
struct Handler<'a> {
    options: ReadOptions,
    cv: Option<&'a ControlledVocabulary>,
    exp: TargetedExperiment,
    warnings: Vec<String>,
    open_tags: Vec<String>,
    publication: Publication,
    contact: Contact,
    instrument: Instrument,
    prediction: Prediction,
    software: Software,
    protein: Protein,
    rt: RetentionTime,
    peptide: Peptide,
    compound: Compound,
    transition: ReactionMonitoringTransition,
    target: IncludeExcludeTarget,
    validation: CVTermList,
    interpretation: Interpretation,
    product: TraMLProduct,
    configuration: Configuration,
    sourcefile: SourceFile,
    sequence_text: Option<String>,
}

impl Handler<'_> {
    /// A message the source logs without losing data.
    fn warn(&mut self, message: String) {
        self.warnings.push(message);
    }

    /// A message the source logs while dropping data: an error by default,
    /// a warning in source-compatible mode.
    fn drop_data(&mut self, message: String) -> Result<()> {
        if self.options.source_compatible {
            self.warnings.push(message);
            Ok(())
        } else {
            Err(Error::Unsupported(format!(
                "{message} (the source drops this data; read with \
                 ReadOptions::source_compatible to do the same)"
            )))
        }
    }

    fn parents(&self) -> (String, String) {
        let n = self.open_tags.len();
        let parent = if n > 1 {
            self.open_tags[n - 2].clone()
        } else {
            String::new()
        };
        let grandparent = if n > 2 {
            self.open_tags[n - 3].clone()
        } else {
            String::new()
        };
        (parent, grandparent)
    }

    fn start(&mut self, tag: &str, attributes: &Attributes) -> Result<()> {
        self.open_tags.push(tag.to_owned());
        if tag == "Sequence" {
            self.sequence_text = Some(String::new());
        }
        if START_IGNORED.contains(&tag) {
            return Ok(());
        }
        let (parent, grandparent) = self.parents();
        match tag {
            "cvParam" => {
                let value = attributes.optional("value").unwrap_or("");
                let unit_accession = attributes.optional("unitAccession").unwrap_or("");
                let unit_name = attributes.optional("unitName").unwrap_or("");
                let unit_cv_ref = attributes.optional("unitCvRef").unwrap_or("");
                let cv_ref = attributes.optional("cvRef").unwrap_or("");
                let mut term = CVTerm::new(
                    attributes.required("accession")?,
                    attributes.required("name")?,
                    cv_ref,
                );
                let mut meta = MetaValue::from(value);
                // Source `CVTerm::hasUnit` tests the unit accession only.
                if !unit_accession.is_empty() {
                    meta = meta.with_unit(Unit::new(unit_accession, unit_name, unit_cv_ref))?;
                }
                term.value = meta;
                self.handle_cv_param(&grandparent, &parent, term)?;
            }
            "userParam" => {
                let xsd_type = attributes.optional("type").unwrap_or("");
                let value = attributes.optional("value").unwrap_or("");
                let name = attributes.required("name")?.to_owned();
                self.handle_user_param(&grandparent, &parent, &name, xsd_type, value)?;
            }
            "cv" => {
                self.exp.add_cv(CV::new(
                    attributes.required("id")?,
                    attributes.required("fullName")?,
                    attributes.required("version")?,
                    attributes.required("URI")?,
                ));
            }
            "Contact" => self.contact.id = attributes.required("id")?.to_owned(),
            "Publication" => self.publication.id = attributes.required("id")?.to_owned(),
            "Instrument" => self.instrument.id = attributes.required("id")?.to_owned(),
            "Software" => {
                self.software.name = attributes.required("id")?.to_owned();
                self.software.version = attributes.required("version")?.to_owned();
            }
            "Protein" => {
                self.protein = Protein::default();
                self.protein.id = attributes.required("id")?.to_owned();
            }
            "Peptide" => {
                self.peptide = Peptide::default();
                self.peptide.id = attributes.required("id")?.to_owned();
                self.peptide.sequence = attributes.required("sequence")?.to_owned();
            }
            "Modification" => {
                // Zero means "no value", as in the source.
                let avg = match attributes.optional("averageMassDelta") {
                    Some(text) => source_double(text)?,
                    None => 0.0,
                };
                let mono = match attributes.optional("monoisotopicMassDelta") {
                    Some(text) => source_double(text)?,
                    None => 0.0,
                };
                // TraML counts locations from 1.
                let location = source_int32(attributes.required("location")?)?
                    .checked_sub(1)
                    .ok_or_else(|| parse_error("Modification location out of range"))?;
                self.peptide.mods.push(PeptideModification {
                    avg_mass_delta: avg,
                    mono_mass_delta: mono,
                    location,
                    ..PeptideModification::default()
                });
            }
            "Compound" => {
                self.compound = Compound::default();
                self.compound.id = attributes.required("id")?.to_owned();
            }
            "Prediction" => {
                self.prediction.software_ref = attributes.required("softwareRef")?.to_owned();
                if let Some(contact) = attributes.optional("contactRef") {
                    self.prediction.contact_ref = contact.to_owned();
                }
            }
            "RetentionTime" => {
                self.rt = RetentionTime::default();
                if let Some(software) = attributes.optional("softwareRef") {
                    self.rt.software_ref = software.to_owned();
                }
            }
            "Transition" => {
                self.transition = ReactionMonitoringTransition::default();
                if let Some(id) = attributes.optional("id") {
                    self.transition.set_name(id);
                }
                if let Some(peptide) = attributes.optional("peptideRef") {
                    self.transition.set_peptide_ref(peptide);
                }
                if let Some(compound) = attributes.optional("compoundRef") {
                    self.transition.set_compound_ref(compound);
                }
            }
            "Interpretation" => {
                if let Some(primary) = attributes.optional("primary") {
                    self.interpretation
                        .cv_terms
                        .metadata
                        .insert("primary".into(), MetaValue::from(primary));
                }
            }
            "Configuration" => {
                self.configuration.instrument_ref =
                    attributes.required("instrumentRef")?.to_owned();
                if let Some(contact) = attributes.optional("contactRef") {
                    self.configuration.contact_ref = contact.to_owned();
                }
            }
            "SourceFile" => {
                self.sourcefile.native_id_type = attributes.required("id")?.to_owned();
                self.sourcefile.name = attributes.required("name")?.to_owned();
                self.sourcefile.path = attributes.required("location")?.to_owned();
            }
            "ProteinRef" => {
                let reference = attributes.required("ref")?.to_owned();
                self.peptide.protein_refs.push(reference);
            }
            "Target" => {
                self.target = IncludeExcludeTarget::default();
                if let Some(id) = attributes.optional("id") {
                    self.target.set_name(id);
                }
                if let Some(peptide) = attributes.optional("peptideRef") {
                    self.target.set_peptide_ref(peptide);
                }
                if let Some(compound) = attributes.optional("compoundRef") {
                    self.target.set_compound_ref(compound);
                }
            }
            _ => {
                self.drop_data(format!("TraMLHandler: unknown tag opening: '{tag}'"))?;
            }
        }
        Ok(())
    }

    fn characters(&mut self, text: &str) {
        if let Some(sequence) = self.sequence_text.as_mut() {
            sequence.push_str(text);
        }
    }

    fn end(&mut self, tag: &str) -> Result<()> {
        let (parent, grandparent) = self.parents();
        self.open_tags.pop();
        if tag == "Sequence" {
            // Source `onCharacters` assigns the protein sequence; an empty
            // element produces no characters and leaves it unchanged.
            if let Some(text) = self.sequence_text.take() {
                if !text.is_empty() {
                    self.protein.sequence = text;
                }
            }
        }
        if END_IGNORED.contains(&tag) {
            return Ok(());
        }
        match tag {
            "Contact" => self.exp.add_contact(std::mem::take(&mut self.contact)),
            "Instrument" => self.exp.add_instrument(std::mem::take(&mut self.instrument)),
            "Publication" => self
                .exp
                .add_publication(std::mem::take(&mut self.publication)),
            "Software" => self.exp.add_software(std::mem::take(&mut self.software)),
            // The source does not reset the pending protein here; the next
            // `<Protein>` start does.
            "Protein" => self.exp.add_protein(self.protein.clone()),
            "RetentionTime" => {
                let rt = std::mem::take(&mut self.rt);
                if grandparent == "Peptide" {
                    self.peptide.rts.push(rt);
                } else if grandparent == "Compound" {
                    self.compound.rts.push(rt);
                } else if parent == "Target" {
                    self.target.set_retention_time(rt);
                } else if parent == "Transition" {
                    self.transition.set_retention_time(rt);
                } else {
                    self.drop_data(format!(
                        "TraMLHandler: tag 'RetentionTime' not allowed at parent tag '{parent}', ignoring!"
                    ))?;
                }
            }
            "Peptide" => self.exp.add_peptide(std::mem::take(&mut self.peptide)),
            "Compound" => self.exp.add_compound(std::mem::take(&mut self.compound)),
            "Transition" => self
                .exp
                .add_transition(std::mem::take(&mut self.transition)),
            "Product" => self
                .transition
                .set_product(std::mem::take(&mut self.product)),
            "IntermediateProduct" => self
                .transition
                .add_intermediate_product(std::mem::take(&mut self.product)),
            "Interpretation" => self
                .product
                .add_interpretation(std::mem::take(&mut self.interpretation)),
            "Prediction" => self
                .transition
                .set_prediction(std::mem::take(&mut self.prediction)),
            "Configuration" => {
                let configuration = std::mem::take(&mut self.configuration);
                if grandparent == "IntermediateProduct" || grandparent == "Product" {
                    self.product.add_configuration(configuration);
                } else if grandparent == "Target" {
                    self.target.add_configuration(configuration);
                } else {
                    self.drop_data(format!(
                        "TraMLHandler: tag 'Configuration' not allowed at parent tag '{parent}', ignoring!"
                    ))?;
                }
            }
            "ValidationStatus" => {
                let validation = std::mem::take(&mut self.validation);
                self.configuration.validations.push(validation);
            }
            "SourceFile" => self
                .exp
                .add_source_file(std::mem::take(&mut self.sourcefile)),
            "Target" => {
                if parent == "TargetIncludeList" {
                    self.exp
                        .add_include_target(std::mem::take(&mut self.target));
                } else if parent == "TargetExcludeList" {
                    self.exp
                        .add_exclude_target(std::mem::take(&mut self.target));
                } else {
                    self.target = IncludeExcludeTarget::default();
                    self.drop_data(format!(
                        "TraMLHandler: tag 'Target' not allowed at parent tag '{parent}', ignoring!"
                    ))?;
                }
            }
            _ => {
                self.drop_data(format!("TraMLHandler: unknown tag closing: '{tag}'"))?;
            }
        }
        Ok(())
    }

    /// The vocabulary checks at the top of the source `handleCVParam_`.
    ///
    /// Returns `false` when the source returns early and the term is not
    /// stored. As in the source, the value-type checks are only reached for
    /// an *obsolete* term (they sit inside the second `if (term.obsolete)`),
    /// and only accessions the loaded PSI-MS vocabulary defines are checked.
    fn check_cv_term(&mut self, parent: &str, term: &CVTerm) -> Result<bool> {
        let Some(cv) = self.cv else {
            return Ok(true);
        };
        // The source loads only psi-ms.obo into this handler's vocabulary.
        if !term.accession.starts_with("MS:") {
            return Ok(true);
        }
        let Ok(definition) = cv.get_term(&term.accession) else {
            return Ok(true);
        };
        let accession = &term.accession;
        let name = definition.name.clone();
        if definition.obsolete {
            self.warn(format!(
                "Obsolete CV term '{accession} - {name}' used in tag '{parent}'."
            ));
        }
        let ws = [' ', '\t', '\n', '\r'];
        let parsed = term.name.trim_matches(ws);
        let correct = definition.name.trim_matches(ws);
        if parsed != correct {
            self.warn(format!(
                "Name of CV term not correct: '{} - {parsed}' should be '{correct}'",
                definition.id
            ));
        }
        if !definition.obsolete {
            return Ok(true);
        }
        self.warn(format!(
            "Obsolete CV term '{accession} - {name}' used in tag '{parent}'."
        ));
        let value = data_value_text(&term.value);
        let xref = definition.xref_type;
        if !value.is_empty() {
            if xref == XRefType::None {
                if !accession.starts_with("PATO:") {
                    self.warn(format!(
                        "The CV term '{accession} - {name}' used in tag '{parent}' must not have a value. The value is '{value}'."
                    ));
                }
                return Ok(true);
            }
            let failure = match xref {
                XRefType::String => None,
                XRefType::Integer
                | XRefType::NegativeInteger
                | XRefType::PositiveInteger
                | XRefType::NonNegativeInteger
                | XRefType::NonPositiveInteger => source_int32(&value)
                    .err()
                    .map(|_| "must have an integer value"),
                XRefType::Decimal => source_double(&value)
                    .err()
                    .map(|_| "must have a floating-point value"),
                XRefType::Date => crate::data_structures::DateTime::parse(&value)
                    .err()
                    .map(|_| "must be a valid date"),
                other => {
                    self.warn(format!(
                        "The CV term '{accession} - {name}' used in tag '{parent}' has the unknown value type '{}'.",
                        other.name()
                    ));
                    None
                }
            };
            if let Some(what) = failure {
                self.drop_data(format!(
                    "The CV term '{accession} - {name}' used in tag '{parent}' {what}. The value is '{value}'."
                ))?;
                return Ok(false);
            }
        } else if xref != XRefType::None && xref != XRefType::String {
            self.drop_data(format!(
                "The CV term '{accession} - {name}' used in tag '{parent}' should have a numerical value. The value is '{value}'."
            ))?;
            return Ok(false);
        }
        Ok(true)
    }

    fn handle_cv_param(&mut self, grandparent: &str, parent: &str, term: CVTerm) -> Result<()> {
        if !self.check_cv_term(parent, &term)? {
            return Ok(());
        }
        let value = data_value_text(&term.value);
        let accession = term.accession.as_str();
        match parent {
            "Software" => self.software.cv_terms.add(term)?,
            "Publication" => self.publication.cv_terms.add(term)?,
            "Instrument" => self.instrument.cv_terms.add(term)?,
            "Contact" => self.contact.cv_terms.add(term)?,
            "RetentionTime" => {
                let unit = term.value.unit().map(|u| u.accession().to_owned());
                match unit.as_deref() {
                    Some("UO:0000010") => self.rt.retention_time_unit = RTUnit::Second,
                    Some("UO:0000031") => self.rt.retention_time_unit = RTUnit::Minute,
                    _ => {
                        // Do not overwrite a unit set by an earlier term.
                        if self.rt.retention_time_unit == RTUnit::SizeOfRTUnit {
                            self.rt.retention_time_unit = RTUnit::Unknown;
                        }
                    }
                }
                match accession {
                    "MS:1000895" => {
                        self.rt.set_rt(source_double(&value)?);
                        self.rt.retention_time_type = RTType::Local;
                    }
                    "MS:1000896" => {
                        self.rt.set_rt(source_double(&value)?);
                        self.rt.retention_time_type = RTType::Normalized;
                    }
                    "MS:1000897" => {
                        self.rt.set_rt(source_double(&value)?);
                        self.rt.retention_time_type = RTType::Predicted;
                    }
                    "MS:1000902" => {
                        if !value.is_empty() {
                            self.rt.set_rt(source_double(&value)?);
                        }
                        self.rt.retention_time_type = RTType::Hpins;
                    }
                    "MS:1002005" => {
                        if !value.is_empty() {
                            self.rt.set_rt(source_double(&value)?);
                        }
                        self.rt.retention_time_type = RTType::Irt;
                    }
                    _ => {
                        self.warn(format!(
                            "The CV term '{accession}' - '{}' used in tag '{parent}' is currently not supported!",
                            term.name
                        ));
                        self.rt.cv_terms.add(term)?;
                    }
                }
            }
            "Evidence" => self.peptide.evidence.add(term)?,
            "Peptide" => match accession {
                "MS:1000041" => self.peptide.set_charge_state(source_int32(&value)?),
                "MS:1000893" => self.peptide.set_peptide_group_label(value),
                "MS:1002476" => self.peptide.set_drift_time(source_double(&value)?),
                _ => self.peptide.cv_terms.add(term)?,
            },
            "Modification" => {
                let modification = self.peptide.mods.last_mut().ok_or_else(|| {
                    parse_error("cvParam inside a Modification that was not opened")
                })?;
                let is_unimod = accession.len() > 7
                    && accession
                        .get(..7)
                        .is_some_and(|p| p.eq_ignore_ascii_case("unimod:"));
                if is_unimod {
                    modification.unimod_id = source_int32(&accession[7..])?;
                } else {
                    modification.cv_terms.add(term)?;
                }
            }
            "Compound" => match accession {
                "MS:1001117" => self.compound.theoretical_mass = source_double(&value)?,
                "MS:1000866" => self.compound.molecular_formula = value,
                "MS:1000868" => self.compound.smiles_string = value,
                "MS:1000041" => self.compound.set_charge_state(source_int32(&value)?),
                // The source sets the drift time of the pending *peptide*
                // here, which the next `<Peptide>` start discards; the
                // compound's own drift time is set instead.
                "MS:1002476" => self.compound.set_drift_time(source_double(&value)?),
                _ => self.compound.cv_terms.add(term)?,
            },
            "Protein" => self.protein.cv_terms.add(term)?,
            "Configuration" => self.configuration.cv_terms.add(term)?,
            "Prediction" => self.prediction.cv_terms.add(term)?,
            "Interpretation" => {
                let ion = match accession {
                    "MS:1000903" => {
                        self.interpretation.ordinal = self.small_count(&value, "ordinal")?;
                        return Ok(());
                    }
                    "MS:1000926" => {
                        self.interpretation.rank = self.small_count(&value, "rank")?;
                        return Ok(());
                    }
                    "MS:1001229" => PeptideFragmentType::AIon,
                    "MS:1001224" => PeptideFragmentType::BIon,
                    "MS:1001231" => PeptideFragmentType::CIon,
                    "MS:1001228" => PeptideFragmentType::XIon,
                    "MS:1001220" => PeptideFragmentType::YIon,
                    "MS:1001230" => PeptideFragmentType::ZIon,
                    "MS:1001523" => PeptideFragmentType::Precursor,
                    "MS:1001222" => PeptideFragmentType::BIonMinusH2O,
                    "MS:1001223" => PeptideFragmentType::YIonMinusH2O,
                    "MS:1001232" => PeptideFragmentType::BIonMinusNH3,
                    "MS:1001233" => PeptideFragmentType::YIonMinusNH3,
                    "MS:1001240" => PeptideFragmentType::NonIdentified,
                    _ => {
                        self.interpretation.cv_terms.add(term)?;
                        return Ok(());
                    }
                };
                self.interpretation.iontype = ion;
            }
            "ValidationStatus" => self.validation.add(term)?,
            "TargetList" => self.exp.add_target_cv_term(term)?,
            "Target" => self.target.cv_terms.add(term)?,
            "Precursor" => {
                if grandparent == "Transition" {
                    if accession == "MS:1000827" {
                        self.transition.set_precursor_mz(source_double(&value)?);
                    } else {
                        self.transition.add_precursor_cv_term(term)?;
                    }
                } else if grandparent == "Target" {
                    self.target.add_precursor_cv_term(term)?;
                } else {
                    // The source ignores this term without a message.
                    self.drop_data(format!(
                        "The CV term '{accession}' used in a Precursor under '{grandparent}' is ignored"
                    ))?;
                }
            }
            "IntermediateProduct" | "Product" => match accession {
                "MS:1000041" => {
                    // The source converts through `toDouble` into an `int`.
                    let charge = source_double(&value)?;
                    if !charge.is_finite() || charge < f64::from(i32::MIN) || charge > f64::from(i32::MAX) {
                        return Err(parse_error(format!(
                            "product charge '{value}' is outside the int range"
                        )));
                    }
                    self.product.set_charge_state(charge.trunc() as i32);
                }
                "MS:1000827" => self.product.set_mz(source_double(&value)?),
                _ => self.product.cv_terms.add(term)?,
            },
            "SourceFile" => self.sourcefile.cv_terms.add(term)?,
            "Transition" => match accession {
                "MS:1002007" => self
                    .transition
                    .set_decoy_transition_type(DecoyTransitionType::Target),
                "MS:1002008" => self
                    .transition
                    .set_decoy_transition_type(DecoyTransitionType::Decoy),
                "MS:1001226" | "MS:1000905" => self
                    .transition
                    .set_library_intensity(source_double(&value)?),
                _ => self.transition.cv_terms.add(term)?,
            },
            _ => {
                let name = term.name.clone();
                self.drop_data(format!(
                    "The CV term '{accession}' - '{name}' used in tag '{parent}' could not be handled, ignoring it!"
                ))?;
            }
        }
        Ok(())
    }

    /// An interpretation ordinal or rank: the source converts with `toInt32`
    /// and stores into an `unsigned char`, wrapping modulo 256. That wrap is
    /// refused by default and reproduced in source-compatible mode.
    fn small_count(&mut self, value: &str, what: &str) -> Result<u8> {
        let number = source_int32(value)?;
        match u8::try_from(number) {
            Ok(small) => Ok(small),
            Err(_) => {
                self.drop_data(format!(
                    "Interpretation {what} {number} does not fit the stored 8-bit range and wraps"
                ))?;
                Ok(number as u8)
            }
        }
    }

    fn handle_user_param(
        &mut self,
        grandparent: &str,
        parent: &str,
        name: &str,
        xsd_type: &str,
        value: &str,
    ) -> Result<()> {
        let data = from_xsd_string(xsd_type, value)?;
        let set = |meta: &mut MetaInfo| {
            meta.insert(name.to_owned(), data.clone());
        };
        match parent {
            "Software" => set(&mut self.software.cv_terms.metadata),
            "Publication" => set(&mut self.publication.cv_terms.metadata),
            "Instrument" => set(&mut self.instrument.cv_terms.metadata),
            "Contact" => set(&mut self.contact.cv_terms.metadata),
            "RetentionTime" => set(&mut self.rt.cv_terms.metadata),
            "Evidence" => set(&mut self.peptide.evidence.metadata),
            "Peptide" => set(&mut self.peptide.cv_terms.metadata),
            "Modification" => {
                let modification = self.peptide.mods.last_mut().ok_or_else(|| {
                    parse_error("userParam inside a Modification that was not opened")
                })?;
                set(&mut modification.cv_terms.metadata);
            }
            "Compound" => set(&mut self.compound.cv_terms.metadata),
            "Protein" => set(&mut self.protein.cv_terms.metadata),
            "Configuration" => set(&mut self.configuration.cv_terms.metadata),
            "Prediction" => set(&mut self.prediction.cv_terms.metadata),
            "Interpretation" => set(&mut self.interpretation.cv_terms.metadata),
            "ValidationStatus" => set(&mut self.validation.metadata),
            "TargetList" => self.exp.set_target_meta_value(name, data),
            "Target" => set(&mut self.target.cv_terms.metadata),
            "Precursor" => {
                // Both land on the owning object, not on its precursor list.
                if grandparent == "Transition" {
                    set(&mut self.transition.cv_terms.metadata);
                } else if grandparent == "Target" {
                    set(&mut self.target.cv_terms.metadata);
                } else {
                    self.drop_data(format!(
                        "userParam '{name}' in a Precursor under '{grandparent}' is ignored"
                    ))?;
                }
            }
            // A product's userParams land on the *transition*, as in the source.
            "Product" => set(&mut self.transition.cv_terms.metadata),
            "SourceFile" => set(&mut self.sourcefile.cv_terms.metadata),
            "Transition" => match name {
                "detecting_transition" => self
                    .transition
                    .set_detecting_transition(value == "true" || value == "1"),
                "identifying_transition" => self
                    .transition
                    .set_identifying_transition(value == "true" || value == "1"),
                "quantifying_transition" => self
                    .transition
                    .set_quantifying_transition(value == "true" || value == "1"),
                _ => set(&mut self.transition.cv_terms.metadata),
            },
            _ => {
                self.drop_data(format!(
                    "Unhandled userParam '{name}' in tag '{parent}'."
                ))?;
            }
        }
        Ok(())
    }
}

/// Attribute-value normalisation as a validating parser applies it: tab,
/// line feed and carriage return become spaces before entities are resolved.
fn attributes_of(start: &BytesStart<'_>) -> Result<Attributes> {
    let mut out: Vec<(String, String)> = Vec::new();
    for attribute in start.attributes() {
        let attribute = attribute.map_err(|e| parse_error(e.to_string()))?;
        let key = std::str::from_utf8(attribute.key.as_ref())
            .map_err(|e| parse_error(e.to_string()))?
            .to_owned();
        let raw =
            std::str::from_utf8(&attribute.value).map_err(|e| parse_error(e.to_string()))?;
        let normalized = raw.replace(['\r', '\n', '\t'], " ");
        let value = quick_xml::escape::unescape(&normalized)
            .map_err(|e| parse_error(e.to_string()))?
            .into_owned();
        if out.iter().any(|(existing, _)| *existing == key) {
            return Err(parse_error(format!("duplicate XML attribute {key}")));
        }
        out.push((key, value));
    }
    Ok(Attributes(out))
}

/// Read a TraML document into a [`TargetedExperiment`]; the source
/// `TraMLHandler` in load mode, driven by `XMLFile::parse_`.
///
/// `cv`, when given, is the vocabulary the source loads into the handler
/// (`psi-ms.obo`) for its term checks; `None` skips those checks and their
/// warnings, which never changes the experiment in default mode.
///
/// # Errors
///
/// [`Error::Parse`] for XML that is not well formed, a document type
/// declaration, a missing required attribute, a number that does not
/// convert, or a ceiling of `options` that is reached;
/// [`Error::Unsupported`] for data the source drops (see the module
/// documentation), unless [`ReadOptions::source_compatible`] is set;
/// [`Error::InvalidValue`] for a CV term that [`CVTermList::add`] refuses
/// (an empty accession or one with whitespace), which the source stores;
/// [`Error::Io`] from `input`.
pub fn read_document(
    input: impl BufRead,
    options: &ReadOptions,
    cv: Option<&ControlledVocabulary>,
) -> Result<TraMLLoad> {
    if options.max_elements == 0 || options.max_depth == 0 {
        return Err(Error::InvalidValue("TraML read limits must be nonzero".into()));
    }
    let limit = options.max_input_bytes;
    let mut reader = Reader::from_reader(input.take(limit.saturating_add(1)));
    reader.config_mut().expand_empty_elements = true;
    reader.config_mut().enable_all_checks(true);
    let mut handler = Handler {
        options: *options,
        cv,
        ..Handler::default()
    };
    let mut buffer = Vec::new();
    let mut elements = 0usize;
    let mut seen_root = false;
    loop {
        let event = reader
            .read_event_into(&mut buffer)
            .map_err(|e| parse_error(e.to_string()))?;
        if reader.buffer_position() > limit {
            return Err(parse_error(format!(
                "TraML document exceeds the {limit}-byte input limit"
            )));
        }
        match event {
            Event::Start(start) => {
                if handler.open_tags.is_empty() && seen_root {
                    return Err(parse_error("multiple XML roots"));
                }
                seen_root = true;
                elements += 1;
                if elements > options.max_elements || handler.open_tags.len() >= options.max_depth
                {
                    return Err(parse_error("TraML element or depth limit exceeded"));
                }
                let name = std::str::from_utf8(start.name().as_ref())
                    .map_err(|e| parse_error(e.to_string()))?
                    .to_owned();
                if handler.open_tags.is_empty() && name != "TraML" {
                    return Err(parse_error(format!("expected a TraML root, found '{name}'")));
                }
                let attributes = attributes_of(&start)?;
                handler.start(&name, &attributes)?;
            }
            Event::End(end) => {
                let name = std::str::from_utf8(end.name().as_ref())
                    .map_err(|e| parse_error(e.to_string()))?
                    .to_owned();
                handler.end(&name)?;
            }
            Event::Text(text) => {
                let decoded = text.xml_content().map_err(|e| parse_error(e.to_string()))?;
                if handler.open_tags.is_empty() {
                    if !decoded.trim_matches([' ', '\t', '\n', '\r']).is_empty() {
                        return Err(parse_error("text outside the TraML root"));
                    }
                } else {
                    handler.characters(&decoded);
                }
            }
            Event::GeneralRef(reference) => {
                let name = reference.decode().map_err(|e| parse_error(e.to_string()))?;
                let mut utf8 = [0u8; 4];
                let text: &str = match reference
                    .resolve_char_ref()
                    .map_err(|e| parse_error(e.to_string()))?
                {
                    Some(c) => c.encode_utf8(&mut utf8),
                    None => quick_xml::escape::resolve_xml_entity(&name)
                        .ok_or_else(|| Error::Unsupported("external XML entity".into()))?,
                };
                handler.characters(text);
            }
            Event::CData(data) => {
                let text = data.decode().map_err(|e| parse_error(e.to_string()))?;
                handler.characters(&text);
            }
            Event::Comment(_) | Event::PI(_) | Event::Decl(_) => {}
            Event::DocType(_) => {
                return Err(Error::Unsupported("DTD in a TraML document".into()));
            }
            Event::Empty(_) => {
                return Err(parse_error("unexpanded empty element"));
            }
            Event::Eof => break,
        }
        buffer.clear();
    }
    if !handler.open_tags.is_empty() {
        return Err(parse_error("truncated TraML document"));
    }
    if !seen_root {
        return Err(parse_error("missing TraML root"));
    }
    Ok(TraMLLoad {
        experiment: handler.exp,
        warnings: handler.warnings,
    })
}

/// The document writer: `TraMLHandler::writeTo` and its `write*_` helpers.
struct Writer<'w> {
    out: &'w mut dyn Write,
    db: &'w ModificationsDB,
}

fn g(value: f64) -> String {
    ostream_g(value, 15)
}

impl Writer<'_> {
    fn line(&mut self, text: &str) -> Result<()> {
        self.out.write_all(text.as_bytes())?;
        self.out.write_all(b"\n")?;
        Ok(())
    }

    /// `writeUserParam_`: one `userParam` per metadata entry.
    ///
    /// The source lists keys in `MetaInfoRegistry` order, which depends on
    /// the order names were first registered in the process; this port's
    /// metadata is ordered by name.
    fn user_params(&mut self, meta: &MetaInfo, indent: usize) -> Result<()> {
        for (key, value) in meta {
            let xsd = match value.data() {
                MetaValueData::Integer(_) => "xsd:integer",
                MetaValueData::Float(_) => "xsd:double",
                _ => "xsd:string",
            };
            let line = format!(
                "{}<userParam name=\"{}\" type=\"{xsd}\" value=\"{}\"/>",
                " ".repeat(2 * indent),
                write_xml_escape(key),
                write_xml_escape(&data_value_text(value))
            );
            self.line(&line)?;
        }
        Ok(())
    }

    /// `writeCVParams_` / `writeCVList_`: every term, by accession.
    fn cv_params(&mut self, terms: &CVTermList, indent: usize) -> Result<()> {
        for list in terms.terms().values() {
            for term in list {
                let mut line = format!(
                    "{}<cvParam cvRef=\"{}\" accession=\"{}\" name=\"{}\"",
                    " ".repeat(2 * indent),
                    attribute_minimal_escape(&term.cv_ref),
                    attribute_minimal_escape(&term.accession),
                    attribute_minimal_escape(&term.name)
                );
                let value = data_value_text(&term.value);
                if term.has_value() && !value.is_empty() {
                    line.push_str(&format!(" value=\"{}\"", attribute_minimal_escape(&value)));
                }
                if let Some(unit) = term.value.unit() {
                    line.push_str(&format!(
                        " unitCvRef=\"{}\" unitAccession=\"{}\" unitName=\"{}\"",
                        attribute_minimal_escape(unit.cv_ref()),
                        attribute_minimal_escape(unit.accession()),
                        attribute_minimal_escape(unit.name())
                    ));
                }
                line.push_str("/>");
                self.line(&line)?;
            }
        }
        Ok(())
    }

    /// CV terms, then userParams: the pair every element body starts with.
    fn params(&mut self, terms: &CVTermList, indent: usize) -> Result<()> {
        self.cv_params(terms, indent)?;
        self.user_params(&terms.metadata, indent)
    }

    fn document(&mut self, exp: &TargetedExperiment) -> Result<()> {
        self.line(r#"<?xml version="1.0" encoding="UTF-8"?>"#)?;
        self.line(r#"<TraML version="1.0.0" xmlns="http://psi.hupo.org/ms/traml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://psi.hupo.org/ms/traml TraML1.0.0.xsd">"#)?;
        self.line("  <cvList>")?;
        if exp.cvs().is_empty() {
            self.line(r#"    <cv id="MS" fullName="Proteomics Standards Initiative Mass Spectrometry Ontology" version="unknown" URI="http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"/>"#)?;
            self.line(r#"    <cv id="UO" fullName="Unit Ontology" version="unknown" URI="http://obo.cvs.sourceforge.net/obo/obo/ontology/phenotype/unit.obo"/>"#)?;
        } else {
            for cv in exp.cvs() {
                let line = format!(
                    "    <cv id=\"{}\" fullName=\"{}\" version=\"{}\" URI=\"{}\"/>",
                    write_xml_escape(&cv.id),
                    write_xml_escape(&cv.fullname),
                    write_xml_escape(&cv.version),
                    write_xml_escape(&cv.uri)
                );
                self.line(&line)?;
            }
        }
        self.line("  </cvList>")?;

        if !exp.source_files().is_empty() {
            self.line("  <SourceFileList>")?;
            for file in exp.source_files() {
                let line = format!(
                    "    <SourceFile id=\"{}\" name=\"{}\" location=\"{}\">",
                    write_xml_escape(&file.native_id_type),
                    write_xml_escape(&file.name),
                    write_xml_escape(&file.path)
                );
                self.line(&line)?;
                self.params(&file.cv_terms, 3)?;
                self.line("    </SourceFile>")?;
            }
            self.line("  </SourceFileList>")?;
        }
        self.simple_list(
            "ContactList",
            "Contact",
            exp.contacts().iter().map(|c| (&c.id, &c.cv_terms)),
        )?;
        self.simple_list(
            "PublicationList",
            "Publication",
            exp.publications().iter().map(|c| (&c.id, &c.cv_terms)),
        )?;
        self.simple_list(
            "InstrumentList",
            "Instrument",
            exp.instruments().iter().map(|c| (&c.id, &c.cv_terms)),
        )?;
        if !exp.software().is_empty() {
            self.line("  <SoftwareList>")?;
            for software in exp.software() {
                let line = format!(
                    "    <Software id=\"{}\" version=\"{}\">",
                    write_xml_escape(&software.name),
                    write_xml_escape(&software.version)
                );
                self.line(&line)?;
                self.params(&software.cv_terms, 3)?;
                self.line("    </Software>")?;
            }
            self.line("  </SoftwareList>")?;
        }
        if !exp.proteins().is_empty() {
            self.line("  <ProteinList>")?;
            for protein in exp.proteins() {
                self.line(&format!(
                    "    <Protein id=\"{}\">",
                    write_xml_escape(&protein.id)
                ))?;
                self.params(&protein.cv_terms, 3)?;
                // The source writes the sequence unescaped.
                self.line(&format!(
                    "      <Sequence>{}</Sequence>",
                    text_minimal_escape(&protein.sequence)
                ))?;
                self.line("    </Protein>")?;
            }
            self.line("  </ProteinList>")?;
        }
        if !exp.compounds().is_empty() || !exp.peptides().is_empty() {
            self.line("  <CompoundList>")?;
            for peptide in exp.peptides() {
                self.peptide(peptide)?;
            }
            for compound in exp.compounds() {
                self.compound(compound)?;
            }
            self.line("  </CompoundList>")?;
        }
        if !exp.transitions().is_empty() {
            self.line("  <TransitionList>")?;
            for transition in exp.transitions() {
                self.transition(transition)?;
            }
            self.line("  </TransitionList>")?;
        }
        if !exp.include_targets().is_empty() || !exp.exclude_targets().is_empty() {
            self.line("  <TargetList>")?;
            self.params(exp.target_cv_terms(), 2)?;
            if !exp.include_targets().is_empty() {
                self.line("    <TargetIncludeList>")?;
                for target in exp.include_targets() {
                    self.target(target)?;
                }
                self.line("    </TargetIncludeList>")?;
            }
            if !exp.exclude_targets().is_empty() {
                self.line("    <TargetExcludeList>")?;
                for target in exp.exclude_targets() {
                    self.target(target)?;
                }
                self.line("    </TargetExcludeList>")?;
            }
            self.line("  </TargetList>")?;
        }
        self.line("</TraML>")
    }

    fn simple_list<'a>(
        &mut self,
        list: &str,
        element: &str,
        items: impl ExactSizeIterator<Item = (&'a String, &'a CVTermList)>,
    ) -> Result<()> {
        if items.len() == 0 {
            return Ok(());
        }
        self.line(&format!("  <{list}>"))?;
        for (id, terms) in items {
            self.line(&format!("    <{element} id=\"{}\">", write_xml_escape(id)))?;
            self.params(terms, 3)?;
            self.line(&format!("    </{element}>"))?;
        }
        self.line(&format!("  </{list}>"))
    }

    fn peptide(&mut self, peptide: &Peptide) -> Result<()> {
        self.line(&format!(
            "    <Peptide id=\"{}\" sequence=\"{}\">",
            write_xml_escape(&peptide.id),
            attribute_minimal_escape(&peptide.sequence)
        ))?;
        if let Some(charge) = peptide.charge_state() {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{charge}"/>"#
            ))?;
        }
        if !peptide.peptide_group_label().is_empty() {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1000893" name="peptide group label" value="{}"/>"#,
                attribute_minimal_escape(peptide.peptide_group_label())
            ))?;
        }
        if peptide.drift_time() >= 0.0 {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1002476" name="ion mobility drift time" value="{}" unitAccession="UO:0000028" unitName="millisecond" unitCvRef="UO" />"#,
                g(peptide.drift_time())
            ))?;
        }
        self.params(&peptide.cv_terms, 3)?;
        for reference in &peptide.protein_refs {
            self.line(&format!(
                "      <ProteinRef ref=\"{}\"/>",
                write_xml_escape(reference)
            ))?;
        }
        for modification in &peptide.mods {
            let mut line = format!(
                "      <Modification location=\"{}\"",
                i64::from(modification.location) + 1
            );
            if modification.mono_mass_delta != 0.0 {
                line.push_str(&format!(
                    " monoisotopicMassDelta=\"{}\"",
                    g(modification.mono_mass_delta)
                ));
            }
            if modification.avg_mass_delta != 0.0 {
                line.push_str(&format!(
                    " averageMassDelta=\"{}\"",
                    g(modification.avg_mass_delta)
                ));
            }
            line.push('>');
            self.line(&line)?;
            if modification.unimod_id != -1 {
                let name = self.unimod_name(peptide, modification)?;
                self.line(&format!(
                    r#"        <cvParam cvRef="UNIMOD" accession="UNIMOD:{}" name="{}"/>"#,
                    modification.unimod_id,
                    attribute_minimal_escape(&name)
                ))?;
            }
            self.params(&modification.cv_terms, 4)?;
            self.line("      </Modification>")?;
        }
        if !peptide.rts.is_empty() {
            self.line("      <RetentionTimeList>")?;
            for rt in &peptide.rts {
                self.retention_time(rt)?;
            }
            self.line("      </RetentionTimeList>")?;
        }
        if !peptide.evidence.is_empty() {
            self.line("      <Evidence>")?;
            self.params(&peptide.evidence, 4)?;
            self.line("      </Evidence>")?;
        }
        self.line("    </Peptide>")
    }

    /// The UniMod short name the source looks up with
    /// `ModificationsDB::getModification("UniMod:<id>", residue, term)`.
    fn unimod_name(&self, peptide: &Peptide, modification: &PeptideModification) -> Result<String> {
        let sequence: Vec<char> = peptide.sequence.chars().collect();
        let (term, residue) = if modification.location < 0 {
            (TermSpecificity::NTerm, sequence.first().copied())
        } else if usize::try_from(modification.location).unwrap_or(usize::MAX) >= sequence.len()
        {
            (TermSpecificity::CTerm, sequence.last().copied())
        } else {
            (
                TermSpecificity::Anywhere,
                usize::try_from(modification.location)
                    .ok()
                    .and_then(|i| sequence.get(i).copied()),
            )
        };
        let name = format!("UniMod:{}", modification.unimod_id);
        // The source picks the first match and only warns when there are more.
        self.db
            .find(&name, residue, Some(term))
            .first()
            .map(|m| m.name().to_owned())
            .ok_or_else(|| {
                Error::InvalidValue(format!(
                    "Retrieving the modification failed. It is not available for the residue '{}' and term specificity '{}'. ({name})",
                    residue.map(String::from).unwrap_or_default(),
                    term.name()
                ))
            })
    }

    fn compound(&mut self, compound: &Compound) -> Result<()> {
        self.line(&format!(
            "    <Compound id=\"{}\">",
            write_xml_escape(&compound.id)
        ))?;
        if let Some(charge) = compound.charge_state() {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{charge}"/>"#
            ))?;
        }
        if compound.theoretical_mass > 0.0 {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1001117" name="theoretical mass" value="{}" unitCvRef="UO" unitAccession="UO:0000221" unitName="dalton"/>"#,
                g(compound.theoretical_mass)
            ))?;
        }
        if !compound.molecular_formula.is_empty() {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1000866" name="molecular formula" value="{}"/>"#,
                attribute_minimal_escape(&compound.molecular_formula)
            ))?;
        }
        if !compound.smiles_string.is_empty() {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1000868" name="SMILES string" value="{}"/>"#,
                attribute_minimal_escape(&compound.smiles_string)
            ))?;
        }
        self.params(&compound.cv_terms, 3)?;
        if !compound.rts.is_empty() {
            self.line("      <RetentionTimeList>")?;
            for rt in &compound.rts {
                self.retention_time(rt)?;
            }
            self.line("      </RetentionTimeList>")?;
        }
        self.line("    </Compound>")
    }

    fn transition(&mut self, transition: &ReactionMonitoringTransition) -> Result<()> {
        let mut line = format!(
            "    <Transition id=\"{}\"",
            write_xml_escape(transition.name())
        );
        if !transition.peptide_ref().is_empty() {
            line.push_str(&format!(
                " peptideRef=\"{}\"",
                write_xml_escape(transition.peptide_ref())
            ));
        }
        if !transition.compound_ref().is_empty() {
            line.push_str(&format!(
                " compoundRef=\"{}\"",
                write_xml_escape(transition.compound_ref())
            ));
        }
        line.push('>');
        self.line(&line)?;
        self.line("      <Precursor>")?;
        self.line(&format!(
            r#"        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
            crate::param::value::format_float(transition.precursor_mz(), true)
        ))?;
        if let Some(precursor) = transition.precursor_cv_term_list() {
            self.params(precursor, 4)?;
        }
        self.line("      </Precursor>")?;
        for product in transition.intermediate_products() {
            self.line("      <IntermediateProduct>")?;
            self.product(product)?;
            self.line("      </IntermediateProduct>")?;
        }
        self.line("      <Product>")?;
        self.product(transition.product())?;
        self.line("      </Product>")?;
        let rt = transition.retention_time();
        if !rt.cv_terms.is_empty() {
            self.retention_time(rt)?;
        }
        if let Ok(prediction) = transition.prediction() {
            let mut line = format!(
                "      <Prediction softwareRef=\"{}\"",
                write_xml_escape(&prediction.software_ref)
            );
            if !prediction.contact_ref.is_empty() {
                line.push_str(&format!(
                    " contactRef=\"{}\"",
                    write_xml_escape(&prediction.contact_ref)
                ));
            }
            line.push('>');
            self.line(&line)?;
            self.params(&prediction.cv_terms, 4)?;
            self.line("      </Prediction>")?;
        }
        self.cv_params(&transition.cv_terms, 3)?;
        if transition.library_intensity() > -100.0 {
            self.line(&format!(
                r#"      <cvParam cvRef="MS" accession="MS:1001226" name="product ion intensity" value="{}"/>"#,
                g(transition.library_intensity())
            ))?;
        }
        match transition.decoy_transition_type() {
            DecoyTransitionType::Target => self.line(
                r#"      <cvParam cvRef="MS" accession="MS:1002007" name="target SRM transition"/>"#,
            )?,
            DecoyTransitionType::Decoy => self.line(
                r#"      <cvParam cvRef="MS" accession="MS:1002008" name="decoy SRM transition"/>"#,
            )?,
            DecoyTransitionType::Unknown => {}
        }
        // Only non-default flags are written; the default is detecting,
        // not identifying, quantifying.
        if !transition.is_detecting_transition() {
            self.line(r#"      <userParam name="detecting_transition" type="xsd:boolean" value="false"/>"#)?;
        }
        if transition.is_identifying_transition() {
            self.line(r#"      <userParam name="identifying_transition" type="xsd:boolean" value="true"/>"#)?;
        }
        if !transition.is_quantifying_transition() {
            self.line(r#"      <userParam name="quantifying_transition" type="xsd:boolean" value="false"/>"#)?;
        }
        self.user_params(&transition.cv_terms.metadata, 3)?;
        self.line("    </Transition>")
    }

    fn retention_time(&mut self, rt: &RetentionTime) -> Result<()> {
        let mut line = String::from("        <RetentionTime");
        if !rt.software_ref.is_empty() {
            line.push_str(&format!(
                " softwareRef=\"{}\"",
                write_xml_escape(&rt.software_ref)
            ));
        }
        // The source opens the value's cvParam only when a value is set and
        // closes it with the unit suffix (or `/>`) unconditionally, so an
        // unset value leaves a stray `/>` or unit text in the element: the
        // document is not well formed. That is refused before anything of
        // the element is written.
        let Some(value) = rt.rt() else {
            return Err(Error::InvalidValue(
                "a retention time without a value cannot be written as TraML".into(),
            ));
        };
        line.push('>');
        self.line(&line)?;
        let (accession, name) = match rt.retention_time_type {
            RTType::Normalized => ("MS:1000896", "normalized retention time"),
            RTType::Predicted => ("MS:1000897", "predicted retention time"),
            RTType::Hpins => ("MS:1000902", "H-PINS retention time normalization standard"),
            RTType::Irt => ("MS:1002005", "iRT retention time normalization standard"),
            _ => ("MS:1000895", "local retention time"),
        };
        let mut line = format!(
            r#"          <cvParam cvRef="MS" accession="{accession}" name="{name}" value="{}""#,
            g(value)
        );
        match rt.retention_time_unit {
            RTUnit::Second => {
                line.push_str(r#" unitCvRef="UO" unitAccession="UO:0000010" unitName="second"/>"#)
            }
            RTUnit::Minute => {
                line.push_str(r#" unitCvRef="UO" unitAccession="UO:0000031" unitName="minute"/>"#)
            }
            _ => line.push_str("/>"),
        }
        self.line(&line)?;
        self.params(&rt.cv_terms, 5)?;
        self.line("        </RetentionTime>")
    }

    fn target(&mut self, target: &IncludeExcludeTarget) -> Result<()> {
        let mut line = format!("      <Target id=\"{}\"", write_xml_escape(target.name()));
        if !target.peptide_ref().is_empty() {
            line.push_str(&format!(
                " peptideRef=\"{}\"",
                write_xml_escape(target.peptide_ref())
            ));
        }
        if !target.compound_ref().is_empty() {
            line.push_str(&format!(
                " compoundRef=\"{}\"",
                write_xml_escape(target.compound_ref())
            ));
        }
        line.push('>');
        self.line(&line)?;
        self.line("        <Precursor>")?;
        self.params(target.precursor_cv_term_list(), 5)?;
        self.line("        </Precursor>")?;
        let rt = target.retention_time();
        if !rt.cv_terms.is_empty() {
            self.retention_time(rt)?;
        }
        if !target.configurations().is_empty() {
            self.line("        <ConfigurationList>")?;
            for configuration in target.configurations() {
                self.configuration(configuration)?;
            }
            self.line("        </ConfigurationList>")?;
        }
        // The source leaves the target's own CV terms and userParams
        // unwritten ("TODO"); so does this port.
        self.line("      </Target>")
    }

    fn product(&mut self, product: &TraMLProduct) -> Result<()> {
        if let Some(charge) = product.charge_state() {
            self.line(&format!(
                r#"        <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="{charge}"/>"#
            ))?;
        }
        if product.mz() > 0.0 {
            self.line(&format!(
                r#"        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="{}" unitCvRef="MS" unitAccession="MS:1000040" unitName="m/z"/>"#,
                g(product.mz())
            ))?;
        }
        self.params(&product.cv_terms, 4)?;
        if !product.interpretation_list().is_empty() {
            self.line("        <InterpretationList>")?;
            for interpretation in product.interpretation_list() {
                self.line("          <Interpretation>")?;
                if interpretation.ordinal > 0 {
                    self.line(&format!(
                        r#"            <cvParam cvRef="MS" accession="MS:1000903" name="product ion series ordinal" value="{}"/>"#,
                        interpretation.ordinal
                    ))?;
                }
                if interpretation.rank > 0 {
                    self.line(&format!(
                        r#"            <cvParam cvRef="MS" accession="MS:1000926" name="product interpretation rank" value="{}"/>"#,
                        interpretation.rank
                    ))?;
                }
                if let Some((accession, name)) = ion_type_term(interpretation.iontype) {
                    self.line(&format!(
                        r#"            <cvParam cvRef="MS" accession="{accession}" name="{name}"/>"#
                    ))?;
                }
                self.params(&interpretation.cv_terms, 6)?;
                self.line("          </Interpretation>")?;
            }
            self.line("        </InterpretationList>")?;
        }
        if !product.configuration_list().is_empty() {
            self.line("        <ConfigurationList>")?;
            for configuration in product.configuration_list() {
                self.configuration(configuration)?;
            }
            self.line("        </ConfigurationList>")?;
        }
        Ok(())
    }

    fn configuration(&mut self, configuration: &Configuration) -> Result<()> {
        let mut line = format!(
            "          <Configuration instrumentRef=\"{}\"",
            write_xml_escape(&configuration.instrument_ref)
        );
        if !configuration.contact_ref.is_empty() {
            line.push_str(&format!(
                " contactRef=\"{}\"",
                write_xml_escape(&configuration.contact_ref)
            ));
        }
        line.push('>');
        self.line(&line)?;
        self.params(&configuration.cv_terms, 6)?;
        for validation in &configuration.validations {
            // Only a validation with CV terms is written, as in the source.
            if !validation.is_empty() {
                self.line("            <ValidationStatus>")?;
                self.params(validation, 7)?;
                self.line("            </ValidationStatus>")?;
            }
        }
        self.line("          </Configuration>")
    }
}

/// Text content escaping for the unescaped `<Sequence>`: `&` and `<` only.
fn text_minimal_escape(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;")
}

/// The PSI-MS term `writeProduct_` writes for an ion type, or `None` for
/// the types it writes nothing for (`Unannotated`, the terminal and full
/// types, and the `Zp1`/`Zp2` ions, for which the source also logs
/// "not supported").
pub fn ion_type_term(ion: PeptideFragmentType) -> Option<(&'static str, &'static str)> {
    Some(match ion {
        PeptideFragmentType::AIon => ("MS:1001229", "frag: a ion"),
        PeptideFragmentType::BIon => ("MS:1001224", "frag: b ion"),
        PeptideFragmentType::CIon => ("MS:1001231", "frag: c ion"),
        PeptideFragmentType::XIon => ("MS:1001228", "frag: x ion"),
        PeptideFragmentType::YIon => ("MS:1001220", "frag: y ion"),
        PeptideFragmentType::ZIon => ("MS:1001230", "frag: z ion"),
        PeptideFragmentType::Precursor => ("MS:1001523", "frag: precursor ion"),
        PeptideFragmentType::BIonMinusH2O => ("MS:1001222", "frag: b ion - H2O"),
        PeptideFragmentType::YIonMinusH2O => ("MS:1001223", "frag: y ion - H2O"),
        PeptideFragmentType::BIonMinusNH3 => ("MS:1001232", "frag: b ion - NH3"),
        PeptideFragmentType::YIonMinusNH3 => ("MS:1001233", "frag: y ion - NH3"),
        PeptideFragmentType::NonIdentified => ("MS:1001240", "non-identified ion"),
        PeptideFragmentType::Unannotated
        | PeptideFragmentType::Zp1Ion
        | PeptideFragmentType::Zp2Ion
        | PeptideFragmentType::Full
        | PeptideFragmentType::Internal
        | PeptideFragmentType::NTerminal
        | PeptideFragmentType::CTerminal => return None,
    })
}

/// Write `exp` as a TraML 1.0.0 document; the source `TraMLHandler` in store
/// mode (`writeTo`), with numbers as `XMLFile::save_` configures the stream
/// (15 significant digits) and the precursor m/z at full precision.
///
/// UniMod modification names are looked up in `db`, as the source looks
/// them up in `ModificationsDB::getInstance()`.
///
/// # Errors
///
/// [`Error::InvalidValue`] for a modification whose UniMod id `db` does not
/// know for its residue and terminus (the source throws
/// `Exception::InvalidValue`), and for a retention time that has CV terms
/// but no value, for which the source writes a malformed element;
/// [`Error::Io`] from `out`.
pub fn write_document(
    out: &mut dyn Write,
    exp: &TargetedExperiment,
    db: &ModificationsDB,
) -> Result<()> {
    Writer { out, db }.document(exp)
}
