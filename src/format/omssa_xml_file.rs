// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Reader for OMSSA XML result files.
//!
//! Ports `FORMAT/OMSSAXMLFile.h`/`.cpp`: loads the documents that implement
//! the OMSSA XML schema (`<MSResponse>` holding one `MSHitSet` per spectrum,
//! each with `MSHits` peptide hits). The source is a Xerces SAX handler; this
//! port drives the same callbacks through the crate's Xerces-faithful SAX
//! driver (`format::sax_xml`), with the limits of
//! [`crate::format::search_engine_xml::XmlReadLimits`].
//! OMSSA reports modifications by its own internal numbers, which the bundled
//! `resources/search_engines/OMSSA_modification_mapping` (the source's
//! `share/OpenMS/CHEMISTRY/OMSSA_modification_mapping`) maps to UniMod names.
//! See `docs/SEARCH_ENGINE_FILES_SUPPORT.md`.

use crate::chemistry::AASequence;
use crate::chemistry::{
    ModificationDefinitionsSet, ModificationsDB, ResidueModification, TermSpecificity,
};
use crate::data_structures::DateTime;
use crate::data_structures::string_utils::{split, to_double, to_int32, trim};
use crate::format::sax_xml::{SaxAttributes as Attributes, SaxEncoding, SaxHandler, drive};
use crate::format::search_engine_text::parse_permissive_sequence;
use crate::format::search_engine_xml::XmlReadLimits;
use crate::identification::{
    FlankingResidue, PeptideEvidence, PeptideHit, PeptideIdentification, ProteinHit,
    ProteinIdentification,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

/// Score type and identifier prefix the reader records.
pub const SCORE_TYPE: &str = "OMSSA";

/// The bundled OMSSA modification mapping, the source's
/// `share/OpenMS/CHEMISTRY/OMSSA_modification_mapping` at the core pin.
pub const BUNDLED_MODIFICATION_MAPPING: &str =
    include_str!("../../resources/search_engines/OMSSA_modification_mapping");

/// The first OMSSA number given to a user-defined variable modification, as
/// the OMSSA adapter numbers `-mux` user modifications.
pub const FIRST_USER_MODIFICATION: u32 = 119;

/// What [`OMSSAXMLFile::load`] reads; the source's two defaulted flags.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OMSSALoadOptions {
    /// Collect a protein hit for every accession the peptide hits reference.
    /// When false, protein accessions are not even stored on the peptide
    /// evidences. Default true.
    pub load_proteins: bool,
    /// Report identifications of spectra without any hit. Default true.
    pub load_empty_hits: bool,
}

impl Default for OMSSALoadOptions {
    fn default() -> Self {
        Self {
            load_proteins: true,
            load_empty_hits: true,
        }
    }
}

/// Everything one OMSSA XML document yields.
#[derive(Clone, Debug)]
pub struct OMSSAXmlResult {
    /// One protein hit per referenced accession, in accession order (none
    /// unless proteins are loaded), score type `OMSSA`, lower is better.
    pub protein_identification: ProteinIdentification,
    /// One identification per `MSHitSet`, in document order, hits sorted by
    /// E-value (lower is better).
    pub peptide_identifications: Vec<PeptideIdentification>,
    /// The source's `warning()` reports, in document order: an ambiguous or
    /// unmapped OMSSA modification number.
    pub warnings: Vec<String>,
}

/// Reader for OMSSA XML files.
///
/// Construction reads the modification mapping, as the source constructor
/// does; [`OMSSAXMLFile::set_modification_definitions`] then registers the
/// search's own modifications.
#[derive(Clone, Debug)]
pub struct OMSSAXMLFile {
    /// Resource ceilings for one document.
    pub limits: XmlReadLimits,
    /// Date of the run, which names the identifier `OMSSA_<date time>` and is
    /// stored as the protein identification's date. `None` takes the local
    /// clock, as the source's `DateTime::now()` does.
    pub run_date: Option<DateTime>,
    mods_map: BTreeMap<u32, Vec<ResidueModification>>,
    mods_to_num: BTreeMap<String, u32>,
    mod_def_set: ModificationDefinitionsSet,
}

impl OMSSAXMLFile {
    /// A reader with the bundled modification mapping.
    ///
    /// # Errors
    ///
    /// As [`OMSSAXMLFile::with_mapping`]; the bundled mapping resolves in the
    /// crate's modification database, so an error here is a packaging defect.
    pub fn new() -> Result<Self> {
        Self::with_mapping(BUNDLED_MODIFICATION_MAPPING)
    }

    /// A reader with a caller's modification mapping text (the source's
    /// private `readMappingFile_`).
    ///
    /// Every line is split at commas; lines that are empty or start with `#`
    /// are skipped. Field one is the OMSSA number, field two a description,
    /// and every further non-blank field a modification name looked up in
    /// [`ModificationsDB::global`].
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] when a line has a single field ("Invalid mapping file
    /// line") or its number does not convert; the lookup's error for an
    /// unknown modification name.
    pub fn with_mapping(mapping: &str) -> Result<Self> {
        let mut mods_map: BTreeMap<u32, Vec<ResidueModification>> = BTreeMap::new();
        let mut mods_to_num = BTreeMap::new();
        let database = ModificationsDB::global();
        for (index, line) in mapping.lines().enumerate() {
            let fields = split(line, ',');
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parse_error = |message: String| Error::Parse {
                line: index + 1,
                message,
            };
            let number = to_int32(trim(fields[0])).map_err(parse_error)?;
            if fields.len() < 2 {
                return Err(parse_error(format!("Invalid mapping file line: '{line}'")));
            }
            let mut mods = Vec::new();
            for name in fields.iter().skip(2).map(|field| trim(field)) {
                if !name.is_empty() {
                    let modification = database.get_modification(name, None, None)?.clone();
                    mods_to_num.insert(modification.full_id().to_owned(), number as u32);
                    mods.push(modification);
                }
            }
            mods_map.insert(number as u32, mods);
        }
        Ok(Self {
            limits: XmlReadLimits::default(),
            run_date: None,
            mods_map,
            mods_to_num,
            mod_def_set: ModificationDefinitionsSet::default(),
        })
    }

    /// Set the modifications of the search (source
    /// `setModificationDefinitionsSet`).
    ///
    /// Fixed modifications are applied to every matching residue of every
    /// peptide read afterwards. Each variable modification that the mapping
    /// does not know gets the next free OMSSA number from 119 on, in name
    /// order. As in the source, numbering restarts at 119 on every call while
    /// earlier registrations stay, so a second call with new modifications
    /// adds a second candidate to the same numbers, and the first one wins
    /// with a warning.
    ///
    /// # Errors
    ///
    /// The lookup's error when a variable modification is not in
    /// [`ModificationsDB::global`]; the reader is then unchanged.
    pub fn set_modification_definitions(
        &mut self,
        modifications: &ModificationDefinitionsSet,
    ) -> Result<()> {
        let database = ModificationsDB::global();
        let mut mods_map = self.mods_map.clone();
        let mut mods_to_num = self.mods_to_num.clone();
        let mut number = FIRST_USER_MODIFICATION;
        for name in modifications.variable_names() {
            if let std::collections::btree_map::Entry::Vacant(slot) = mods_to_num.entry(name) {
                let modification = database.get_modification(slot.key(), None, None)?.clone();
                mods_map.entry(number).or_default().push(modification);
                slot.insert(number);
                number = number.wrapping_add(1);
            }
        }
        self.mods_map = mods_map;
        self.mods_to_num = mods_to_num;
        self.mod_def_set = modifications.clone();
        Ok(())
    }

    /// The UniMod candidates of an OMSSA modification number, if mapped.
    pub fn mapped_modifications(&self, number: u32) -> Option<&[ResidueModification]> {
        self.mods_map.get(&number).map(Vec::as_slice)
    }

    /// Load an OMSSA XML file.
    ///
    /// # Errors
    ///
    /// [`Error::Io`] when the file cannot be read; see
    /// [`OMSSAXMLFile::load_reader`] for the rest.
    pub fn load(
        &self,
        path: impl AsRef<Path>,
        options: OMSSALoadOptions,
    ) -> Result<OMSSAXmlResult> {
        self.load_reader(BufReader::new(File::open(path)?), options)
    }

    /// Load OMSSA XML from a reader.
    ///
    /// Every `MSHits` becomes a [`PeptideHit`] with the E-value
    /// (`MSHits_evalue`) as score, the charge of `MSHits_charge`, the
    /// sequence of `MSHits_pepstring` (with the fixed modifications applied
    /// and each `MSModHit` resolved through the mapping) and one
    /// [`PeptideEvidence`] per `MSPepHit`; `MSHits_pepstart` and
    /// `MSHits_pepstop` give the first evidence its flanking residues.
    /// `MSHitSet_ids_E` values of the form `<m/z>_<RT>[_...]` set the
    /// identification's m/z and RT; a value that does not convert is ignored,
    /// as in the source.
    ///
    /// A peptide string that does not parse becomes an empty sequence, as in
    /// the source.
    ///
    /// # Errors
    ///
    /// [`Error::Parse`] for XML that is not well formed, an E-value, charge,
    /// modification site or modification number that does not convert (the
    /// source's `ConversionError`), a modification that cannot be placed, an
    /// invalid flanking residue or a non-finite E-value. [`Error::Unsupported`]
    /// for a DTD or a non-UTF-8 encoding, [`Error::InvalidValue`] when a limit
    /// in [`OMSSAXMLFile::limits`] is exceeded. Nothing is returned on error.
    pub fn load_reader(
        &self,
        reader: impl Read,
        options: OMSSALoadOptions,
    ) -> Result<OMSSAXmlResult> {
        let limits = self.limits.sax()?;
        let mut handler = Handler {
            reader: self,
            options,
            peptide_identifications: Vec::new(),
            actual_peptide_hit: PeptideHit::default(),
            actual_peptide_evidence: PeptideEvidence::default(),
            actual_peptide_evidences: Vec::new(),
            actual_peptide_id: PeptideIdentification::default(),
            tag: String::new(),
            actual_mod_site: 0,
            actual_mod_type: String::new(),
            warnings: Vec::new(),
        };
        drive(
            BufReader::new(reader),
            SaxEncoding::Utf8,
            &limits,
            &mut handler,
        )?;
        let now = self.run_date.unwrap_or_else(DateTime::now);
        let identifier = format!("{SCORE_TYPE}_{}", now.get());
        let mut peptide_identifications = handler.peptide_identifications;
        let mut accessions = BTreeSet::new();
        for identification in &mut peptide_identifications {
            identification.score_type = SCORE_TYPE.into();
            identification.higher_score_better = false;
            identification.identifier = identifier.clone();
            identification.sort()?;
            if options.load_proteins {
                for hit in &identification.hits {
                    accessions.extend(hit.protein_accessions().into_iter().map(str::to_owned));
                }
            }
        }
        let mut protein_identification = ProteinIdentification::default();
        if options.load_proteins {
            protein_identification.hits = accessions
                .into_iter()
                .map(|accession| ProteinHit {
                    accession,
                    ..ProteinHit::default()
                })
                .collect();
            protein_identification.higher_score_better = false;
            protein_identification.score_type = SCORE_TYPE.into();
        }
        protein_identification.date_time = Some(now.get());
        protein_identification.identifier = identifier;
        Ok(OMSSAXmlResult {
            protein_identification,
            peptide_identifications,
            warnings: handler.warnings,
        })
    }
}

struct Handler<'a> {
    reader: &'a OMSSAXMLFile,
    options: OMSSALoadOptions,
    peptide_identifications: Vec<PeptideIdentification>,
    actual_peptide_hit: PeptideHit,
    actual_peptide_evidence: PeptideEvidence,
    actual_peptide_evidences: Vec<PeptideEvidence>,
    actual_peptide_id: PeptideIdentification,
    tag: String,
    actual_mod_site: u32,
    actual_mod_type: String,
    warnings: Vec<String>,
}

fn parse_error(message: impl Into<String>) -> Error {
    Error::Parse {
        line: 0,
        message: message.into(),
    }
}

fn flanking(value: &str) -> Result<Option<FlankingResidue>> {
    match value.chars().next() {
        None => Ok(None),
        Some(code) => FlankingResidue::from_code(code)
            .map(Some)
            .map_err(|_| parse_error(format!("invalid flanking residue '{code}'"))),
    }
}

impl Handler<'_> {
    /// Source `MSHits_pepstring`: a permissive parse (an empty sequence on
    /// failure), then every fixed modification on each residue of its origin.
    fn peptide_string(&mut self, value: &str) -> Result<()> {
        let Ok(mut sequence) = parse_permissive_sequence(trim(value)) else {
            self.actual_peptide_hit.sequence = AASequence::default();
            return Ok(());
        };
        let database = ModificationsDB::global();
        for name in self.reader.mod_def_set.fixed_names() {
            let origin = database.get_modification(&name, None, None)?.origin();
            let positions: Vec<usize> = sequence
                .as_str()
                .chars()
                .enumerate()
                .filter(|(_, residue)| Some(*residue) == origin)
                .map(|(position, _)| position)
                .collect();
            for position in positions {
                sequence.set_modification(position, &name)?;
            }
        }
        self.actual_peptide_hit.sequence = sequence;
        Ok(())
    }

    /// Source `onEndElement` for `MSModHit`.
    fn modification_hit(&mut self) -> Result<()> {
        let number = to_int32(&self.actual_mod_type).map_err(parse_error)? as u32;
        let candidates = self
            .reader
            .mods_map
            .get(&number)
            .filter(|candidates| !candidates.is_empty());
        let Some(candidates) = candidates else {
            self.warnings.push(format!(
                "Cannot find PSI-MOD mapping for mod - ignoring '{}'",
                self.actual_mod_type
            ));
            return Ok(());
        };
        if candidates.len() > 1 {
            self.warnings.push(format!(
                "Cannot determine exact type of modification of position {} in sequence {} using modification {} - using first possibility!",
                self.actual_mod_site, self.actual_peptide_hit.sequence, self.actual_mod_type
            ));
        }
        let modification = &candidates[0];
        let sequence = &mut self.actual_peptide_hit.sequence;
        match modification.term_specificity() {
            TermSpecificity::NTerm => {
                sequence.set_n_terminal_modification(modification.full_id())?
            }
            TermSpecificity::CTerm => {
                sequence.set_c_terminal_modification(modification.full_id())?
            }
            _ => {
                sequence.set_modification(self.actual_mod_site as usize, modification.full_id())?
            }
        }
        Ok(())
    }

    /// Source `MSHitSet_ids_E`: `<m/z>_<RT>[_...]`.
    fn spectrum_ids(&mut self, value: &str) {
        if value.is_empty() || !value.contains('_') {
            return;
        }
        let parts = split(value, '_');
        let Ok(mz) = to_double(parts[0]) else {
            return;
        };
        self.actual_peptide_id.mz = Some(mz);
        if let Ok(rt) = to_double(parts[1]) {
            self.actual_peptide_id.rt = Some(rt);
        }
    }
}

impl SaxHandler for Handler<'_> {
    fn start_element(&mut self, name: &str, _attributes: &Attributes) -> Result<()> {
        self.tag = trim(name).to_owned();
        Ok(())
    }

    fn end_element(&mut self, name: &str) -> Result<()> {
        match trim(name) {
            "MSHits" => {
                self.actual_peptide_hit.evidences =
                    std::mem::take(&mut self.actual_peptide_evidences);
                self.actual_peptide_evidence = PeptideEvidence::default();
                let hit = std::mem::take(&mut self.actual_peptide_hit);
                self.actual_peptide_id.hits.push(hit);
            }
            "MSHitSet" => {
                let identification = std::mem::take(&mut self.actual_peptide_id);
                if !identification.hits.is_empty() || self.options.load_empty_hits {
                    self.peptide_identifications.push(identification);
                }
            }
            "MSModHit" => self.modification_hit()?,
            _ => {}
        }
        self.tag.clear();
        Ok(())
    }

    fn characters(&mut self, text: &str) -> Result<()> {
        if self.tag.is_empty() {
            return Ok(());
        }
        let value = trim(text);
        match self.tag.as_str() {
            "MSPepHit_start"
            | "MSPepHit_stop"
            | "MSPepHit_defline"
            | "MSPepHit_protlength"
            | "MSHits_pvalue"
            | "MSHits_mass"
            | "MSHits_theomass" => {}
            "MSPepHit_accession" => {
                if self.options.load_proteins {
                    self.actual_peptide_evidence.protein_accession = value.to_owned();
                }
            }
            "MSPepHit_oid" => {
                // End of MSPepHit: the evidence is added, and not reset.
                self.actual_peptide_evidences
                    .push(self.actual_peptide_evidence.clone());
            }
            "MSHits_evalue" => {
                let score = to_double(value).map_err(parse_error)?;
                if !score.is_finite() {
                    return Err(parse_error(format!("non-finite E-value {score}")));
                }
                self.actual_peptide_hit.score = score;
            }
            "MSHits_charge" => {
                self.actual_peptide_hit.charge = to_int32(value).map_err(parse_error)?;
            }
            "MSHits_pepstring" => self.peptide_string(value)?,
            "MSHits_pepstart" => {
                if let (Some(code), Some(first)) =
                    (flanking(value)?, self.actual_peptide_evidences.first_mut())
                {
                    first.aa_before = code;
                }
            }
            "MSHits_pepstop" => {
                if let (Some(code), Some(first)) =
                    (flanking(value)?, self.actual_peptide_evidences.first_mut())
                {
                    first.aa_after = code;
                }
            }
            // The modification and spectrum tags keep `tag`, so the white
            // space that follows reaches them again.
            "MSHits_mods" => {
                self.actual_mod_site = 0;
                self.actual_mod_type.clear();
                return Ok(());
            }
            "MSModHit_site" => {
                self.actual_mod_site = to_int32(value).map_err(parse_error)? as u32;
                return Ok(());
            }
            "MSMod" => {
                self.actual_mod_type = value.to_owned();
                return Ok(());
            }
            "MSHitSet_ids_E" => {
                self.spectrum_ids(value);
                return Ok(());
            }
            _ => return Ok(()),
        }
        self.tag.clear();
        Ok(())
    }
}
