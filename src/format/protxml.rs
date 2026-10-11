// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for ProteinProphet's protXML, ported from `FORMAT/ProtXMLFile.h`.
//!
//! Only the parts of the `protein_summary` subtree that carry protein-peptide
//! associations are read; every other element and attribute is ignored. The
//! schema ships with the TPP (`protXML_v6.xsd` in the source's
//! `share/OpenMS/SCHEMAS`); it is not consulted.
//!
//! For protein groups, only the group leader carries its own probability and
//! coverage. Indistinguishable siblings (`<indistinguishable_protein>`) carry
//! only a `protein_name`; on read they inherit the leader's score, so that
//! score-based filtering does not tear groups apart, and they have no coverage.
//!
//! ProteinProphet assigns probability 0 to "unneeded" proteins (subset and
//! subsumable proteins whose peptides are explained by higher-probability ones).
//! Those `<protein>` elements are dropped entirely: no protein hit, no group
//! entry and no peptide evidence is created for them or their children, and a
//! `<protein_group>` left with no accession is not inserted. This keeps
//! downstream tools such as `ProteinQuantifier` from treating subsumable proteins
//! as indistinguishable groups and discarding shared peptides.
//!
//! What a read fills in:
//!
//! - [`ProteinHit`](crate::identification::ProteinHit): accession (`protein_name`), score (the ProteinProphet
//!   probability; siblings inherit the leader's) and coverage
//!   (`percent_coverage`, leader only). No meta values.
//! - [`PeptideHit`](crate::identification::PeptideHit): sequence (`peptide_sequence`, with the modifications of
//!   `<mod_aminoacid_mass>`), score (`nsp_adjusted_probability`), charge, one
//!   [`PeptideEvidence`](crate::identification::PeptideEvidence) per protein of the enclosing indistinguishable group,
//!   and the integer meta values `is_unique` (`is_nondegenerate_evidence` is
//!   `Y`) and `is_contributing` (`is_contributing_evidence` is `Y`).
//! - The score type of both identifications is `ProteinProphet probability`,
//!   higher is better.
//!
//! Writing protXML is not supported, as in the source. The API mapping, the
//! preserved source behaviour and every divergence are listed in
//! `docs/PROTXML_SUPPORT.md`.

use crate::chemistry::{AASequence, ModificationsDB, TermSpecificity};
use crate::concept::unique_id::UniqueIdGenerator;
use crate::data_structures::DateTime;
use crate::data_structures::string_utils::to_double;
use crate::format::pepxml::{find_enzyme, internal_residue_mass};
use crate::format::sax_xml::{
    self, SaxAttributes as Attributes, SaxEncoding, SaxEvent, SaxLimits, xerces_parse_int,
};
use crate::identification::{
    PeptideEvidence, PeptideHit, PeptideIdentification, ProteinGroup, ProteinHit,
    ProteinIdentification,
};
use crate::metadata::MetaValue;
use crate::param::value::format_float;
use crate::{Error, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Score type of both identifications a read produces.
pub const SCORE_TYPE: &str = "ProteinProphet probability";

/// Mass window of the modification look-up in `matchModification_`, in Da.
pub const MODIFICATION_TOLERANCE: f64 = 0.001;

/// Resource ceilings of one read. The defaults accept every upstream fixture
/// and any ProteinProphet output for a single experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReadLimits {
    /// Maximum decoded document size in bytes.
    pub max_bytes: usize,
    /// Maximum XML events, including those of ignored elements.
    pub max_events: usize,
    /// Maximum element nesting depth.
    pub max_depth: usize,
    /// Maximum protein hits (leaders and siblings).
    pub max_proteins: usize,
    /// Maximum peptide hits.
    pub max_peptides: usize,
    /// Maximum retained log lines; later ones are counted in
    /// [`ProtXmlDocument::suppressed_log_lines`].
    pub max_log_lines: usize,
}

impl Default for ReadLimits {
    fn default() -> Self {
        Self {
            max_bytes: 1 << 30,
            max_events: 100_000_000,
            max_depth: 64,
            max_proteins: 10_000_000,
            max_peptides: 10_000_000,
            max_log_lines: 10_000,
        }
    }
}

/// What one read produced: the source's two out-parameters and its log.
#[derive(Clone, Debug, PartialEq)]
pub struct ProtXmlDocument {
    /// Protein hits, protein groups and indistinguishable groups.
    pub protein_identification: ProteinIdentification,
    /// One peptide hit per `<peptide>` element that was not dropped, in
    /// document order; a degenerate peptide appears once per group.
    pub peptide_identification: PeptideIdentification,
    /// The source's warning and error output, in order: missing
    /// `percent_coverage` and `charge`, and modifications that could not be
    /// found, parsed or uniquely identified.
    pub log: Vec<String>,
    /// Log lines beyond [`ReadLimits::max_log_lines`].
    pub suppressed_log_lines: usize,
}

/// The source's `ProtXMLFile`, a stateless reader.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ProtXmlFile;

impl ProtXmlFile {
    /// A reader, the source's constructor.
    pub fn new() -> Self {
        Self
    }

    /// Load a protXML file, as `load(filename, protein_ids, peptide_ids)`.
    ///
    /// Each `<program_details>` element sets the identifier of both
    /// identifications to the next ID of `ids`, in decimal, as the source does
    /// with `UniqueIdGenerator::getUniqueId()`; pass a seeded generator for a
    /// reproducible identifier.
    ///
    /// # Errors
    ///
    /// - I/O errors, and [`Error::InvalidValue`] when a [`ReadLimits`] ceiling
    ///   is exceeded;
    /// - [`Error::Parse`] for a document that is not well-formed XML, and when
    ///   a required attribute is missing or a number does not convert (the
    ///   source's `Exception::ParseError`);
    /// - [`Error::InvalidValue`] when `sample_enzyme` names no enzyme (the
    ///   source's `Exception::ElementNotFound`), when `time` is not a date the
    ///   source's `DateTime::set` accepts (its `Exception::ParseError`), when a
    ///   peptide sequence does not parse, or when a modification the look-up
    ///   found cannot be set on the sequence;
    /// - [`Error::InvalidRange`] when a `<mod_aminoacid_mass>` position is not
    ///   a residue of its peptide (the source's `Exception::IndexOverflow`);
    /// - [`Error::Unsupported`] for a DTD or a reference that needs one.
    ///
    /// Stricter than the source, each an error rather than undefined or
    /// unrepresentable behaviour: an `<indistinguishable_protein>`,
    /// `<peptide>` or `<mod_aminoacid_mass>` with no enclosing protein or
    /// peptide to attach to (the source dereferences an empty vector or a
    /// null or dangling pointer), a charge outside the 32-bit range (the
    /// source's `(int)` cast of `strtol` wraps it), and a probability,
    /// coverage or peptide score that is not finite.
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ProtXmlDocument> {
        self.load_with_limits(path, ids, &ReadLimits::default())
    }

    /// [`Self::load`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn load_with_limits(
        &self,
        path: impl AsRef<Path>,
        ids: &mut UniqueIdGenerator,
        limits: &ReadLimits,
    ) -> Result<ProtXmlDocument> {
        let path = path.as_ref();
        let file = File::open(path)?;
        read_document(
            BufReader::new(file),
            &path.display().to_string(),
            ids,
            limits,
        )
    }

    /// Read protXML from `reader`; see [`Self::load`]. The file name in log
    /// lines is empty.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read(
        &self,
        reader: impl BufRead,
        ids: &mut UniqueIdGenerator,
    ) -> Result<ProtXmlDocument> {
        read_document(reader, "", ids, &ReadLimits::default())
    }

    /// [`Self::read`] with explicit resource ceilings.
    ///
    /// # Errors
    ///
    /// As for [`Self::load`].
    pub fn read_with_limits(
        &self,
        reader: impl BufRead,
        ids: &mut UniqueIdGenerator,
        limits: &ReadLimits,
    ) -> Result<ProtXmlDocument> {
        read_document(reader, "", ids, limits)
    }

    /// Store protXML, as `store(filename, protein_ids, peptide_ids,
    /// document_id)`; not implemented, as in the source.
    ///
    /// # Errors
    ///
    /// Always [`Error::Unsupported`], the source's `Exception::NotImplemented`.
    pub fn store(
        &self,
        _path: impl AsRef<Path>,
        _protein_ids: &ProteinIdentification,
        _peptide_ids: &PeptideIdentification,
        _document_id: &str,
    ) -> Result<()> {
        Err(Error::Unsupported(
            "writing protXML is not implemented, as in the source".into(),
        ))
    }
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

/// `attributeAsString_`: a missing attribute is the source's `fatalError`.
fn required<'a>(attributes: &'a Attributes, name: &str) -> Result<&'a str> {
    attributes.required(name)
}

/// `attributeAsDouble_` / `optionalAttributeAsDouble_`: `StringUtils::toDouble`,
/// whose conversion error is a parse error here, as the source's
/// `ConversionError` ends the parse.
fn double(text: &str, name: &str) -> Result<f64> {
    to_double(text).map_err(|_| {
        parse_error(format!(
            "Could not convert '{text}' in '{name}' to a double"
        ))
    })
}

/// [`double`] for a value the crate stores and requires to be finite.
fn finite(text: &str, name: &str) -> Result<f64> {
    let value = double(text, name)?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::InvalidValue(format!(
            "protXML attribute '{name}' is not finite: {text}"
        )))
    }
}

fn read_document(
    reader: impl BufRead,
    file: &str,
    ids: &mut UniqueIdGenerator,
    limits: &ReadLimits,
) -> Result<ProtXmlDocument> {
    let sax = SaxLimits {
        max_input_bytes: u64::try_from(limits.max_bytes).unwrap_or(u64::MAX),
        max_events: u64::try_from(limits.max_events).unwrap_or(u64::MAX),
        max_depth: limits.max_depth,
        max_attributes: 256,
        // Character data is not read; a run cannot exceed the document.
        max_text_bytes: limits.max_bytes,
    };
    let mut handler = Handler {
        file,
        ids,
        limits: *limits,
        modifications: ModificationsDB::global(),
        protein: ProteinIdentification::default(),
        peptide: PeptideIdentification::default(),
        group: ProteinGroup::default(),
        hit: None,
        skip_protein: false,
        peptides: 0,
        log: Vec::new(),
        suppressed: 0,
    };
    sax_xml::parse_with(reader, SaxEncoding::Declared, &sax, |event| match event {
        SaxEvent::Start(name, attributes) => handler.start(name, attributes),
        SaxEvent::End(name) => handler.end(name),
        SaxEvent::Characters(_) => Ok(()),
    })?;
    Ok(ProtXmlDocument {
        protein_identification: handler.protein,
        peptide_identification: handler.peptide,
        log: handler.log,
        suppressed_log_lines: handler.suppressed,
    })
}

/// The handler state of the source's `ProtXMLFile`, whose members are the
/// `prot_id_`/`pep_id_` targets, `pep_hit_`, `protein_group_` and
/// `skip_protein_`.
struct Handler<'a> {
    file: &'a str,
    ids: &'a mut UniqueIdGenerator,
    limits: ReadLimits,
    modifications: &'static ModificationsDB,
    protein: ProteinIdentification,
    peptide: PeptideIdentification,
    group: ProteinGroup,
    /// `pep_hit_`: the open `<peptide>`; `None` where the source's pointer is
    /// null or dangling.
    hit: Option<PeptideHit>,
    skip_protein: bool,
    peptides: usize,
    log: Vec<String>,
    suppressed: usize,
}

impl Handler<'_> {
    fn log(&mut self, line: String) {
        if self.log.len() < self.limits.max_log_lines {
            self.log.push(line);
        } else {
            self.suppressed += 1;
        }
    }

    /// `XMLHandler::error(LOAD, message)`: logged, not thrown.
    fn error(&mut self, message: &str) {
        let line = format!("Non-fatal error while loading '{}': {message}", self.file);
        self.log(line);
    }

    fn start(&mut self, tag: &str, attributes: &Attributes) -> Result<()> {
        if tag == "protein_summary_header" {
            let database = required(attributes, "reference_database")?;
            let enzyme = required(attributes, "sample_enzyme")?;
            let enzyme = find_enzyme(enzyme).ok_or_else(|| {
                Error::InvalidValue(format!("the enzyme '{enzyme}' was not found"))
            })?;
            self.protein.search_parameters.database = database.to_owned();
            self.protein.search_parameters.digestion_enzyme = enzyme.name().to_owned();
            self.protein.score_type = SCORE_TYPE.into();
            self.protein.higher_score_better = true;
            self.peptide.score_type = SCORE_TYPE.into();
            self.peptide.higher_score_better = true;
        }
        // Not `else if`: the source tests this tag independently.
        if tag == "program_details" {
            let analysis = required(attributes, "analysis")?;
            let time = required(attributes, "time")?;
            let version = required(attributes, "version")?;
            // The source's DateTime::set throws on anything it cannot parse,
            // so its "Cannot parse 'time'" warning is never reached.
            let date = DateTime::parse(time)
                .map_err(|_| Error::InvalidValue(format!("Invalid date time string: '{time}'")))?;
            self.protein.date_time = Some(date.iso_string());
            self.protein.search_engine = analysis.to_owned();
            self.protein.search_engine_version = version.to_owned();
            let identifier = self.ids.get_unique_id().to_string();
            self.protein.identifier = identifier.clone();
            self.peptide.identifier = identifier;
        }

        match tag {
            "protein_group" => {
                self.group = ProteinGroup {
                    probability: finite(required(attributes, "probability")?, "probability")?,
                    ..ProteinGroup::default()
                };
                self.skip_protein = false;
            }
            "protein" => {
                let probability = double(required(attributes, "probability")?, "probability")?;
                // ProteinProphet writes exactly "0.0000" for unneeded
                // proteins, which parses to exactly 0.0.
                if probability == 0.0 {
                    self.skip_protein = true;
                    return Ok(());
                }
                self.skip_protein = false;
                if !probability.is_finite() {
                    return Err(Error::InvalidValue(format!(
                        "protXML protein probability is not finite: {probability}"
                    )));
                }
                let name = required(attributes, "protein_name")?;
                // Open a new indistinguishable group.
                self.protein
                    .indistinguishable_groups
                    .push(ProteinGroup::default());
                self.register_protein(name)?;
                let coverage = match attributes.get("percent_coverage") {
                    Some(text) => Some(finite(text, "percent_coverage")?),
                    None => {
                        self.log("Required attribute 'percent_coverage' missing".into());
                        None
                    }
                };
                let hit = self
                    .protein
                    .hits
                    .last_mut()
                    .expect("register_protein appended a hit");
                hit.coverage = coverage;
                hit.score = probability;
            }
            "indistinguishable_protein" => {
                if self.skip_protein {
                    return Ok(());
                }
                let name = required(attributes, "protein_name")?;
                // The current last protein leads the same group.
                let score = self
                    .protein
                    .hits
                    .last()
                    .map(|hit| hit.score)
                    .ok_or_else(|| {
                        Error::InvalidValue(
                            "protXML indistinguishable_protein precedes every protein".into(),
                        )
                    })?;
                self.register_protein(name)?;
                // The leader's score may not strictly transfer (protein length
                // etc.), but keeping it lets score filters keep groups whole.
                if let Some(hit) = self.protein.hits.last_mut() {
                    hit.score = score;
                }
            }
            "peptide" => {
                if self.skip_protein {
                    return Ok(());
                }
                if self.peptides >= self.limits.max_peptides {
                    return Err(Error::InvalidValue("protXML peptide limit exceeded".into()));
                }
                // A degenerate peptide appears in several groups with
                // different statistics; each instance is its own hit.
                let sequence = AASequence::parse_with_registry(
                    required(attributes, "peptide_sequence")?,
                    self.modifications,
                )?;
                let score = finite(
                    required(attributes, "nsp_adjusted_probability")?,
                    "nsp_adjusted_probability",
                )?;
                let mut hit = PeptideHit {
                    sequence,
                    score,
                    ..PeptideHit::default()
                };
                match attributes.get("charge") {
                    Some(text) => hit.charge = xerces_parse_int(text)?,
                    None => self.log("Required attribute 'charge' missing".into()),
                }
                let group = self
                    .protein
                    .indistinguishable_groups
                    .last()
                    .ok_or_else(|| {
                        Error::InvalidValue("protXML peptide precedes every protein".into())
                    })?;
                hit.evidences = group
                    .accessions
                    .iter()
                    .map(|accession| PeptideEvidence {
                        protein_accession: accession.clone(),
                        ..PeptideEvidence::default()
                    })
                    .collect();
                let unique = required(attributes, "is_nondegenerate_evidence")? == "Y";
                let contributing = required(attributes, "is_contributing_evidence")? == "Y";
                hit.metadata
                    .insert("is_unique".into(), MetaValue::from(i32::from(unique)));
                hit.metadata.insert(
                    "is_contributing".into(),
                    MetaValue::from(i32::from(contributing)),
                );
                self.hit = Some(hit);
            }
            "mod_aminoacid_mass" => {
                if self.skip_protein {
                    return Ok(());
                }
                self.modification(attributes)?;
            }
            _ => {}
        }
        Ok(())
    }

    fn end(&mut self, tag: &str) -> Result<()> {
        match tag {
            "protein_group" => {
                // A group whose proteins were all dropped stays out: consumers
                // assume `accessions[0]` exists.
                if !self.group.accessions.is_empty() {
                    self.protein
                        .protein_groups
                        .push(std::mem::take(&mut self.group));
                }
            }
            "protein" => self.skip_protein = false,
            "peptide" => {
                if self.skip_protein {
                    return Ok(());
                }
                let hit = self.hit.take().ok_or_else(|| {
                    Error::InvalidValue("protXML peptide end without an open peptide".into())
                })?;
                self.peptide.hits.push(hit);
                self.peptides += 1;
            }
            _ => {}
        }
        Ok(())
    }

    /// `registerProtein_`: a new protein hit, appended to the current group
    /// and to the current indistinguishable group.
    fn register_protein(&mut self, name: &str) -> Result<()> {
        if self.protein.hits.len() >= self.limits.max_proteins {
            return Err(Error::InvalidValue("protXML protein limit exceeded".into()));
        }
        let group = self
            .protein
            .indistinguishable_groups
            .last_mut()
            .ok_or_else(|| {
                Error::InvalidValue("protXML protein outside any protein element".into())
            })?;
        group.accessions.push(name.to_owned());
        self.group.accessions.push(name.to_owned());
        self.protein.hits.push(ProteinHit {
            accession: name.to_owned(),
            ..ProteinHit::default()
        });
        Ok(())
    }

    /// The `<mod_aminoacid_mass>` branch, which modifies the last seen peptide.
    fn modification(&mut self, attributes: &Attributes) -> Result<()> {
        let position = xerces_parse_int(required(attributes, "position")?)?;
        let mass = double(required(attributes, "mass")?, "mass")?;
        let mut sequence = self
            .hit
            .as_ref()
            .ok_or_else(|| {
                Error::InvalidValue("protXML mod_aminoacid_mass outside a peptide".into())
            })?
            .sequence
            .clone();
        // `Size position = attributeAsInt_(...)`, then `sequence[position - 1]`:
        // position 0 or a negative one wraps to a huge index, which
        // `AASequence::operator[]` refuses with IndexOverflow like any other
        // index past the end.
        let index = usize::try_from(position)
            .ok()
            .and_then(|position| position.checked_sub(1))
            .filter(|&index| index < sequence.len())
            .ok_or_else(|| {
                Error::InvalidRange(format!(
                    "protXML modification position {position} is outside the {}-residue peptide",
                    sequence.len()
                ))
            })?;
        let origin = char::from(sequence.as_str().as_bytes()[index]);
        let description = self.match_modification(mass, origin)?;
        let position_text = position.to_string();
        match description {
            None => {
                let message = format!(
                    "Cannot find modification '{} {origin}' @{position_text}",
                    format_float(mass, true)
                );
                self.error(&message);
            }
            Some(description) => {
                let parts: Vec<&str> = description.split(' ').collect();
                if parts.len() == 2 {
                    let term = self
                        .modifications
                        .get_modification(&description, None, None)?
                        .term_specificity();
                    if parts[1] == "(C-term)" || term == TermSpecificity::CTerm {
                        sequence.set_c_terminal_modification_with_registry(
                            parts[0],
                            self.modifications,
                        )?;
                    } else if parts[1] == "(N-term)" || term == TermSpecificity::NTerm {
                        sequence.set_n_terminal_modification_with_registry(
                            parts[0],
                            self.modifications,
                        )?;
                    } else {
                        sequence.set_modification_with_registry(
                            index,
                            parts[0],
                            self.modifications,
                        )?;
                    }
                } else {
                    let message =
                        format!("Cannot parse modification '{description}@{position_text}'");
                    self.error(&message);
                }
            }
        }
        if let Some(hit) = self.hit.as_mut() {
            hit.sequence = sequence;
        }
        Ok(())
    }

    /// `matchModification_`: the full ID of the first registry modification
    /// whose monoisotopic difference lies within [`MODIFICATION_TOLERANCE`] of
    /// `mass` less the residue's internal mass, for `origin`, in registry
    /// order. Several candidates are logged and the first is used.
    fn match_modification(&mut self, mass: f64, origin: char) -> Result<Option<String>> {
        let difference = mass - internal_residue_mass(origin)?;
        // A NaN or infinite mass matches nothing in the source's comparisons;
        // the registry search refuses it, so it is answered here.
        if !difference.is_finite() {
            return Ok(None);
        }
        let candidates = self.modifications.search_by_mass(
            difference,
            MODIFICATION_TOLERANCE,
            Some(origin),
            None,
        )?;
        let Some(first) = candidates.first() else {
            return Ok(None);
        };
        let first = first.full_id().to_owned();
        if candidates.len() > 1 {
            let all: Vec<&str> = candidates.iter().map(|m| m.full_id()).collect();
            let message = format!(
                "Modification '{}' is not uniquely defined by the given data. Using '{first}' to represent any of '{}'!",
                format_float(mass, true),
                all.join(", ")
            );
            self.error(&message);
        }
        Ok(Some(first))
    }
}
