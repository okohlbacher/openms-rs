// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for X! Tandem XML result files.
//!
//! Ports `FORMAT/XTandemXMLFile.h`/`.cpp`: loads the documents that
//! implement the X! Tandem output schema (`<bioml>` holding `model` groups of
//! `protein` elements whose `peptide` holds the `domain` peptide hits). The
//! source is a Xerces SAX handler; this port drives the same callbacks through
//! the crate's Xerces-faithful SAX driver (`format::sax_xml`), which splits a
//! text node at every reference as Xerces does, so a `<note>` holding an
//! entity reference keeps only its first run, as in the source (the handler
//! acts on the first `characters` call and then clears its flags). See
//! `docs/SEARCH_ENGINE_FILES_SUPPORT.md` for the API mapping and every
//! divergence.

use crate::chemistry::{
    AASequence, ModificationDefinition, ModificationDefinitionsSet, ModificationMatchOptions,
    ModificationsDB, ResidueModification, TermSpecificity,
};
use crate::data_structures::DateTime;
use crate::data_structures::string_utils::{to_double, to_int32};
use crate::format::sax_xml::{
    SaxAttributes as Attributes, SaxEncoding, SaxHandler, drive, xerces_parse_int,
};
use crate::format::search_engine_xml::XmlReadLimits;
use crate::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification,
};
use crate::metadata::MetaValue;
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Score type and search engine name the reader records.
pub const SCORE_TYPE: &str = "XTandem";

/// The N-terminal modifications X! Tandem applies by default (its parameters
/// "protein, quick pyrolidone" and "protein, quick acetyl"), consulted when
/// an N-terminal mass shift matches none of the caller's definitions.
pub const DEFAULT_N_TERMINAL_MODIFICATIONS: &str =
    "Gln->pyro-Glu (N-term Q),Glu->pyro-Glu (N-term E),Acetyl (N-term)";

/// The mass tolerance, in Da, of every modification lookup.
pub const MODIFICATION_TOLERANCE: f64 = 0.01;

/// Everything one X! Tandem document yields.
#[derive(Clone, Debug)]
pub struct XTandemXmlResult {
    /// Protein hits with the `expect` value as score (lower is better),
    /// sorted ascending, score type and search engine `XTandem`.
    pub protein_identification: ProteinIdentification,
    /// One identification per X! Tandem spectrum id, ordered by that id, each
    /// with its hits sorted by hyperscore (higher is better).
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// The caller's definitions, extended by every X! Tandem default
    /// N-terminal modification the document used (as a variable one).
    pub modification_definitions: ModificationDefinitionsSet,
    /// The source's non-fatal `error()` reports, in document order: a mass
    /// shift that fits no known modification.
    pub warnings: Vec<String>,
}

/// Reader for X! Tandem XML files.
///
/// The document is read as ISO-8859-1 whatever it declares, as the source's
/// `enforceEncoding_("ISO-8859-1")` forces.
#[derive(Clone, Copy, Debug, Default)]
pub struct XTandemXMLFile {
    /// Resource ceilings for one document.
    pub limits: XmlReadLimits,
    /// Date of the run, which names the identifier `XTandem_<date>` and is
    /// stored as the protein identification's date. `None` takes the local
    /// clock, as the source's `DateTime::now()` does; set it to make two
    /// loads of one document equal.
    pub run_date: Option<DateTime>,
}

impl XTandemXMLFile {
    /// A reader with default limits and the local clock as run date.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load an X! Tandem XML file.
    ///
    /// `modifications` are the fixed and variable modifications defined for
    /// the search; a document may extend them with X! Tandem's default
    /// N-terminal modifications, and the extended set is returned in
    /// [`XTandemXmlResult::modification_definitions`].
    ///
    /// The source's `load(filename, protein_identification, id_data,
    /// mod_def_set)` swaps the protein hits into the caller's
    /// `protein_identification` and keeps its other members; the port returns
    /// a fresh one, so a caller that relied on earlier members copies the
    /// hits instead.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read; see
    /// [`XTandemXMLFile::load_reader`] for the rest.
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        modifications: &ModificationDefinitionsSet,
    ) -> Result<XTandemXmlResult> {
        self.load_reader(BufReader::new(File::open(path)?), modifications)
    }

    /// Load X! Tandem XML from a reader; see [`XTandemXMLFile::load`].
    ///
    /// Every `domain` becomes a [`PeptideHit`] of the identification named by
    /// the part of its `id` before the first `.`, with the `hyperscore` as
    /// score and the meta values `nextscore`, `delta`, `mass` (from `mh`),
    /// `E-Value` (from `expect`) and `hyperscore`, plus each present
    /// `<ion>_score` (a double) and `<ion>_ions` (an integer) for the ions a,
    /// b, c, x, y and z. A `domain` repeating the previous sequence of the
    /// same spectrum adds a [`PeptideEvidence`] to that hit instead. `aa`
    /// elements inside a `model` group modify the current hit, `protein`
    /// elements become protein hits (once per `uid`), and the `note` labelled
    /// `Description` names the spectrum reference.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for XML that is not well formed, a missing required
    /// attribute (the source's `fatalError`), a numeric attribute that does not
    /// convert, a sequence that does not parse strictly (the source calls
    /// `AASequence::fromString(seq, false)`), an invalid flanking residue, a
    /// modification whose position lies outside its peptide, or a
    /// non-finite score. Also for an `aa` element outside every `group` or
    /// before any `domain` of its spectrum, where the source logs a
    /// non-fatal error or not at all and then reads an empty container
    /// (undefined behaviour). [`Error::Unsupported`] for a DTD,
    /// [`Error::InvalidValue`] when a limit in [`XTandemXMLFile::limits`] is
    /// exceeded. Nothing is returned on error.
    pub fn load_reader(
        &self,
        reader: impl Read,
        modifications: &ModificationDefinitionsSet,
    ) -> Result<XTandemXmlResult> {
        let limits = self.limits.sax()?;
        let mut handler = Handler::new(modifications.clone())?;
        drive(
            BufReader::new(reader),
            SaxEncoding::Latin1,
            &limits,
            &mut handler,
        )?;
        handler.finish(self.run_date.unwrap_or_else(DateTime::now))
    }
}

/// The source's `GroupType`: the kind of the enclosing `group` element.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GroupType {
    Model,
    Parameters,
    Support,
}

struct Handler {
    mod_def_set: ModificationDefinitionsSet,
    default_nterm_mods: ModificationDefinitionsSet,
    is_protein_note: bool,
    is_spectrum_note: bool,
    skip_protein_acc_update: bool,
    peptide_hits: BTreeMap<u32, Vec<PeptideHit>>,
    protein_hits: Vec<ProteinHit>,
    protein_uids: BTreeSet<u32>,
    current_protein: String,
    current_charge: i32,
    current_id: u32,
    tag: String,
    current_start: u32,
    previous_seq: String,
    spectrum_ids: BTreeMap<u32, String>,
    group_types: Vec<GroupType>,
    warnings: Vec<String>,
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn double(attributes: &Attributes, name: &str) -> Result<f64> {
    to_double(attributes.required(name)?).map_err(parse_error)
}

fn meta_double(value: f64, name: &str) -> Result<MetaValue> {
    MetaValue::try_from(value)
        .map_err(|_| parse_error(format!("non-finite value {value} for '{name}'")))
}

/// The source's `UInt x = attributeAsInt_(...)`: the parsed `Int`, wrapped.
fn unsigned(attributes: &Attributes, name: &str) -> Result<u32> {
    Ok(xerces_parse_int(attributes.required(name)?)? as u32)
}

/// The source's `pe.setStart(position - 1)` on a `UInt`: a one-based
/// position of 0 wraps to -1, the unknown position.
fn zero_based(position: u32) -> Option<usize> {
    let value = position.wrapping_sub(1) as i32;
    usize::try_from(value).ok()
}

fn flanking(code: Option<char>) -> Result<FlankingResidue> {
    match code {
        None => Ok(FlankingResidue::Unknown),
        Some(code) => FlankingResidue::from_code(code)
            .map_err(|_| parse_error(format!("invalid flanking residue '{code}'"))),
    }
}

impl Handler {
    fn new(mod_def_set: ModificationDefinitionsSet) -> Result<Self> {
        Ok(Self {
            mod_def_set,
            default_nterm_mods: ModificationDefinitionsSet::from_comma_separated(
                "",
                DEFAULT_N_TERMINAL_MODIFICATIONS,
            )?,
            is_protein_note: false,
            is_spectrum_note: false,
            skip_protein_acc_update: false,
            peptide_hits: BTreeMap::new(),
            protein_hits: Vec::new(),
            protein_uids: BTreeSet::new(),
            current_protein: String::new(),
            current_charge: 0,
            current_id: 0,
            tag: String::new(),
            current_start: 0,
            previous_seq: String::new(),
            spectrum_ids: BTreeMap::new(),
            group_types: Vec::new(),
            warnings: Vec::new(),
        })
    }

    fn finish(self, now: DateTime) -> Result<XTandemXmlResult> {
        let identifier = format!("XTandem_{}", now.date_string());
        let mut spectrum_ids = self.spectrum_ids;
        let mut peptide_identifications = Vec::with_capacity(self.peptide_hits.len());
        for (id, hits) in self.peptide_hits {
            let mut identification = PeptideIdentification {
                identifier: identifier.clone(),
                hits,
                score_type: SCORE_TYPE.into(),
                higher_score_better: true,
                ..PeptideIdentification::default()
            };
            identification.set_spectrum_reference(spectrum_ids.remove(&id).unwrap_or_default());
            identification.sort()?;
            peptide_identifications.push(identification);
        }
        let mut protein_identification = ProteinIdentification {
            hits: self.protein_hits,
            higher_score_better: false,
            score_type: SCORE_TYPE.into(),
            search_engine: SCORE_TYPE.into(),
            date_time: Some(now.get()),
            identifier,
            ..ProteinIdentification::default()
        };
        // Source `ProteinIdentification::sort` with `ProteinHit::ScoreLess`,
        // which orders by (score, accession): equal E-values fall back to
        // the accession.
        protein_identification.sort()?;
        Ok(XTandemXmlResult {
            protein_identification,
            peptide_identifications,
            modification_definitions: self.mod_def_set,
            warnings: self.warnings,
        })
    }

    fn domain(&mut self, attributes: &Attributes) -> Result<()> {
        let id_string = attributes.required("id")?;
        let prefix = id_string.split('.').next().unwrap_or(id_string);
        let id = to_int32(prefix).map_err(parse_error)? as u32;
        self.current_id = id;

        let aa_before = flanking(attributes.required("pre")?.chars().next_back())?;
        let aa_after = flanking(attributes.required("post")?.chars().next())?;
        self.current_start = unsigned(attributes, "start")?;
        let current_stop = unsigned(attributes, "end")?;
        let evidence = PeptideEvidence {
            protein_accession: self.current_protein.clone(),
            start: zero_based(self.current_start),
            end: zero_based(current_stop),
            aa_before,
            aa_after,
        };

        let seq = attributes.required("seq")?;
        if !self.peptide_hits.contains_key(&id) || seq != self.previous_seq {
            let sequence = AASequence::parse(seq)?;
            let hyperscore = double(attributes, "hyperscore")?;
            let mut hit = PeptideHit {
                sequence,
                charge: self.current_charge,
                ..PeptideHit::default()
            };
            let meta = &mut hit.metadata;
            meta.insert(
                "nextscore".into(),
                meta_double(double(attributes, "nextscore")?, "nextscore")?,
            );
            meta.insert(
                "delta".into(),
                meta_double(double(attributes, "delta")?, "delta")?,
            );
            meta.insert("mass".into(), meta_double(double(attributes, "mh")?, "mh")?);
            meta.insert(
                "E-Value".into(),
                meta_double(double(attributes, "expect")?, "expect")?,
            );
            // A non-finite hyperscore is refused here, as the meta value.
            meta.insert("hyperscore".into(), meta_double(hyperscore, "hyperscore")?);
            hit.score = hyperscore;
            for ion in ['a', 'b', 'c', 'x', 'y', 'z'] {
                let score_name = format!("{ion}_score");
                if let Some(text) = attributes.get(&score_name) {
                    let score = to_double(text).map_err(parse_error)?;
                    let value = meta_double(score, &score_name)?;
                    hit.metadata.insert(score_name, value);
                }
                let count_name = format!("{ion}_ions");
                if let Some(text) = attributes.get(&count_name) {
                    let count = xerces_parse_int(text)? as u32;
                    hit.metadata.insert(count_name, MetaValue::from(count));
                }
            }
            self.peptide_hits.entry(id).or_default().push(hit);
            self.previous_seq = seq.to_owned();
        }
        if let Some(hit) = self
            .peptide_hits
            .get_mut(&id)
            .and_then(|hits| hits.last_mut())
        {
            hit.evidences.push(evidence);
        }
        Ok(())
    }

    /// Candidate definitions for a mass shift, best match first.
    fn find_matches(
        set: &ModificationDefinitionsSet,
        mass_shift: f64,
        aa: &str,
        term: TermSpecificity,
    ) -> Result<Vec<ModificationDefinition>> {
        let options = ModificationMatchOptions {
            residue: aa.to_owned(),
            term_specificity: Some(term),
            tolerance: MODIFICATION_TOLERANCE,
            ..ModificationMatchOptions::default()
        };
        Ok(set
            .find_matches(mass_shift, &options)?
            .into_iter()
            .map(|found| found.definition)
            .collect())
    }

    fn amino_acid(&mut self, attributes: &Attributes) -> Result<()> {
        let Some(group_type) = self.group_types.last().copied() else {
            return Err(parse_error(
                "Found an 'aa' element outside of a 'group'! Please check your input file.",
            ));
        };
        if group_type != GroupType::Model {
            return Ok(());
        }
        let aa = attributes.required("type")?.to_owned();
        let at = xerces_parse_int(attributes.required("at")?)?;
        let mass_shift = double(attributes, "modified")?;
        let mut aa_seq = self
            .peptide_hits
            .get(&self.current_id)
            .and_then(|hits| hits.last())
            .map(|hit| hit.sequence.clone())
            .ok_or_else(|| {
                parse_error("Found an 'aa' element before any 'domain' of its spectrum")
            })?;
        // X! Tandem uses the position in the protein; the source subtracts
        // the `UInt` start with wrapping arithmetic.
        let mod_pos = (at as u32).wrapping_sub(self.current_start) as i32;
        // Source `Int(aa_seq.size() - 1)`, where an empty sequence wraps to -1.
        let last = (aa_seq.len() as i64 - 1) as i32;

        let mut matches = Vec::new();
        if mod_pos <= 0 {
            matches =
                Self::find_matches(&self.mod_def_set, mass_shift, &aa, TermSpecificity::NTerm)?;
        } else if mod_pos >= last {
            matches =
                Self::find_matches(&self.mod_def_set, mass_shift, &aa, TermSpecificity::CTerm)?;
        }
        if matches.is_empty() {
            matches = Self::find_matches(
                &self.mod_def_set,
                mass_shift,
                &aa,
                TermSpecificity::Anywhere,
            )?;
        }
        if matches.is_empty() && mod_pos <= 0 {
            matches = Self::find_matches(
                &self.default_nterm_mods,
                mass_shift,
                &aa,
                TermSpecificity::NTerm,
            )?;
            if let Some(best) = matches.first() {
                let definition =
                    ModificationDefinition::with_options(best.modification_name(), false, 0)?;
                self.mod_def_set.add_modification(definition)?;
            }
        }
        let residue: Option<ResidueModification> = match matches.first() {
            Some(best) => Some(best.modification()?.clone()),
            None => {
                let database = ModificationsDB::global();
                let origin = aa.chars().next();
                let mut found = None;
                if mod_pos <= 0 {
                    found = database.best_by_mass(
                        mass_shift,
                        MODIFICATION_TOLERANCE,
                        origin,
                        Some(TermSpecificity::NTerm),
                    )?;
                } else if mod_pos >= last {
                    found = database.best_by_mass(
                        mass_shift,
                        MODIFICATION_TOLERANCE,
                        origin,
                        Some(TermSpecificity::CTerm),
                    )?;
                }
                if found.is_none() {
                    found = database.best_by_mass(
                        mass_shift,
                        MODIFICATION_TOLERANCE,
                        origin,
                        Some(TermSpecificity::Anywhere),
                    )?;
                }
                found.cloned()
            }
        };
        let Some(modification) = residue else {
            self.warnings.push(format!(
                "No modification found which fits residue '{aa}' with mass '{mass_shift}'!"
            ));
            return Ok(());
        };
        let full_id = modification.full_id();
        match modification.term_specificity() {
            TermSpecificity::NTerm => aa_seq.set_n_terminal_modification(full_id)?,
            TermSpecificity::CTerm => aa_seq.set_c_terminal_modification(full_id)?,
            _ => {
                let index = usize::try_from(mod_pos).map_err(|_| {
                    parse_error(format!(
                        "modification position {mod_pos} lies before the peptide"
                    ))
                })?;
                aa_seq.set_modification(index, full_id)?;
            }
        }
        if let Some(hit) = self
            .peptide_hits
            .get_mut(&self.current_id)
            .and_then(|hits| hits.last_mut())
        {
            hit.sequence = aa_seq;
        }
        Ok(())
    }

    fn group(&mut self, attributes: &Attributes) -> Result<()> {
        match attributes.required("type")? {
            "model" => {
                self.group_types.push(GroupType::Model);
                if let Some(charge) = attributes.get("z") {
                    self.current_charge = to_int32(charge).map_err(parse_error)?;
                }
                self.previous_seq.clear();
            }
            "parameter" => self.group_types.push(GroupType::Parameters),
            _ => self.group_types.push(GroupType::Support),
        }
        Ok(())
    }

    fn protein(&mut self, attributes: &Attributes) -> Result<()> {
        let uid = unsigned(attributes, "uid")?;
        self.current_protein = attributes.required("label")?.to_owned();
        if self.protein_uids.insert(uid) {
            let mut hit = ProteinHit {
                accession: self.current_protein.clone(),
                ..ProteinHit::default()
            };
            if let Some(expect) = attributes.get("expect") {
                let score = to_double(expect).map_err(parse_error)?;
                if !score.is_finite() {
                    return Err(parse_error(format!("non-finite protein expect {score}")));
                }
                hit.score = score;
            }
            self.protein_hits.push(hit);
            self.skip_protein_acc_update = false;
        } else {
            self.skip_protein_acc_update = true;
        }
        Ok(())
    }
}

impl SaxHandler for Handler {
    fn start_element(&mut self, name: &str, attributes: &Attributes) -> Result<()> {
        self.tag = name.to_owned();
        match name {
            "domain" => self.domain(attributes),
            "aa" => self.amino_acid(attributes),
            "group" => self.group(attributes),
            "note" => {
                match attributes.get("label") {
                    Some("description") => self.is_protein_note = true,
                    Some("Description") => self.is_spectrum_note = true,
                    _ => {}
                }
                Ok(())
            }
            "protein" => self.protein(attributes),
            _ => Ok(()),
        }
    }

    fn end_element(&mut self, name: &str) -> Result<()> {
        self.tag = name.to_owned();
        if name == "group" {
            self.group_types.pop();
        }
        Ok(())
    }

    fn characters(&mut self, text: &str) -> Result<()> {
        if self.tag != "note" {
            return Ok(());
        }
        let value = text.trim_matches([' ', '\t', '\n', '\r']);
        if self.is_protein_note {
            self.current_protein = value.to_owned();
            if !self.skip_protein_acc_update {
                if let Some(hit) = self.protein_hits.last_mut() {
                    hit.accession = self.current_protein.clone();
                }
            }
        } else if self.is_spectrum_note {
            self.spectrum_ids.insert(self.current_id, value.to_owned());
        }
        self.is_protein_note = false;
        self.is_spectrum_note = false;
        Ok(())
    }
}
