// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded OBODataProvider translation from OpenMS4-core revision 7c029e8.

use super::{
    EmpiricalFormula, ModificationRecord, ResidueModification, TermSpecificity, full_identifier,
};
use crate::{Error, Result};
use std::collections::BTreeSet;
use std::io::BufRead;

/// Limits include the newline in line/input bytes. Record, alias and registry
/// payload limits apply to the complete registry after an extension.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OboReadOptions {
    pub cross_links_only: bool,
    pub max_input_bytes: usize,
    pub max_line_bytes: usize,
    pub max_terms: usize,
    pub max_records: usize,
    pub max_aliases: usize,
    /// Conservative string, formula, record and index payload bound; allocator
    /// bookkeeping/capacity rounding and transient parser buffers are additional.
    pub max_registry_bytes: usize,
}
impl Default for OboReadOptions {
    fn default() -> Self {
        Self {
            cross_links_only: false,
            max_input_bytes: 16 * 1024 * 1024,
            max_line_bytes: 64 * 1024,
            max_terms: 100_000,
            max_records: 200_000,
            max_aliases: 1_000_000,
            max_registry_bytes: 128 * 1024 * 1024,
        }
    }
}
impl OboReadOptions {
    pub(super) fn validate(&self) -> Result<()> {
        if self.max_line_bytes == 0 {
            return Err(Error::InvalidValue(
                "OBO max_line_bytes must be positive".into(),
            ));
        }
        Ok(())
    }
}

/// Counts for one atomic provider append. Unresolved alias specificity records
/// are omitted, and may share the same OBO accession.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct OboLoadReport {
    pub records_added: usize,
    pub aliases_added: usize,
    pub unresolved_aliases: usize,
}

pub(super) fn checked(current: usize, extra: usize, limit: usize, label: &str) -> Result<usize> {
    current
        .checked_add(extra)
        .filter(|&sum| sum <= limit)
        .ok_or_else(|| Error::InvalidValue(format!("{label} limit exceeded")))
}
fn error(line: usize, message: impl Into<String>) -> Error {
    Error::Parse {
        line,
        message: message.into(),
    }
}

#[derive(Default)]
struct Term {
    record: ModificationRecord,
    accession: String,
    origins: String,
    skip: bool,
    line: usize,
}
struct Parsed {
    records: Vec<(String, ResidueModification)>,
    aliases: usize,
    bytes: usize,
}
impl Parsed {
    fn push_term(&mut self, mut term: Term, options: &OboReadOptions) -> Result<()> {
        if term.skip || term.accession.is_empty() {
            return Ok(());
        }
        let origins: BTreeSet<&str> = term.origins.split(',').collect();
        let mut sites = Vec::new();
        for origin in origins {
            if origin.len() == 1 && !matches!(origin, "B" | "J" | "Z") {
                let residue = origin.as_bytes()[0];
                if !residue.is_ascii_uppercase() {
                    return Err(error(
                        term.line,
                        "OBO origins must be uppercase ASCII residues",
                    ));
                }
                sites.push((Some(residue as char), term.record.term_specificity));
            }
        }
        for (label, peptide, protein) in [
            (
                "ProteinN-term",
                TermSpecificity::NTerm,
                TermSpecificity::ProteinNTerm,
            ),
            (
                "ProteinC-term",
                TermSpecificity::CTerm,
                TermSpecificity::ProteinCTerm,
            ),
        ] {
            if term.origins.contains(label) {
                sites.push((
                    None,
                    if options.cross_links_only {
                        peptide
                    } else {
                        protein
                    },
                ));
            }
        }
        if sites.is_empty() {
            return Ok(());
        }
        term.record.provenance = super::ModificationProvenance::Cv;
        let mut base = ResidueModification::from_record(term.record)
            .map_err(|e| error(term.line, e.to_string()))?;
        for (origin, specificity) in sites {
            let alias = base.record_id.is_some_and(|id| id > 0);
            let wildcard = origin.is_none() || origin == Some('X');
            if !alias
                && wildcard
                && (specificity == TermSpecificity::Anywhere || base.diff_mono_mass == 0.0)
            {
                continue;
            }
            if !alias && base.full_name.is_empty() {
                return Err(error(term.line, "OBO modification requires a name"));
            }
            base.origin = if wildcard && specificity != TermSpecificity::Anywhere {
                None
            } else {
                origin
            };
            base.term = specificity;
            base.full_id = full_identifier(&base.full_name, base.origin, specificity);
            checked(
                self.records.len(),
                1,
                options.max_records,
                "OBO expanded records",
            )?;
            self.aliases = checked(
                self.aliases,
                if alias { 1 } else { base.names().len() },
                options.max_aliases,
                "OBO expanded aliases",
            )?;
            self.bytes = checked(
                self.bytes,
                base.payload_bytes()?,
                options.max_registry_bytes,
                "OBO expanded record bytes",
            )?;
            self.records.push((term.accession.clone(), base.clone()));
        }
        Ok(())
    }
}

pub(super) fn read(
    mut reader: impl BufRead,
    options: &OboReadOptions,
) -> Result<Vec<ResidueModification>> {
    options.validate()?;
    let mut parsed = Parsed {
        records: Vec::new(),
        aliases: 0,
        bytes: 0,
    };
    let mut term = Term::default();
    let mut active = false;
    let mut line = 0;
    let mut terms = 0;
    let mut input_bytes = 0;
    let mut synonym_count = 0;
    while let Some(bytes) = read_line(&mut reader, &mut input_bytes, options)? {
        line += 1;
        let text = std::str::from_utf8(&bytes)
            .map_err(|_| error(line, "OBO input must be UTF-8"))?
            .trim();
        if text
            .chars()
            .any(|c| c.is_control() && !c.is_ascii_whitespace())
        {
            return Err(error(line, "invalid control character in OBO input"));
        }
        if text.is_empty() || text.starts_with('!') {
            continue;
        }
        let compact: String = text.chars().filter(|c| !c.is_whitespace()).collect();
        if compact.starts_with('[') {
            parsed.push_term(std::mem::take(&mut term), options)?;
            active = compact == "[Term]";
            if active {
                terms = checked(terms, 1, options.max_terms, "OBO terms")?;
                term.line = line;
            }
            continue;
        }
        if !active {
            continue;
        }
        let Some((key, value)) = text.split_once(':') else {
            continue;
        };
        let key: String = key.chars().filter(|c| !c.is_whitespace()).collect();
        let value = value.trim();
        match key.as_str() {
            "id" => {
                if value.is_empty() {
                    return Err(error(line, "empty OBO accession"));
                }
                term.accession = value.into();
                term.record.name = value.into();
                term.record.obo_accession = Some(value.into());
            }
            "name" => {
                term.record.full_name = value.into();
                if term.record.name.contains("XLMOD") {
                    term.record.name = value.into();
                }
            }
            "def" => {
                let stripped = value.replace(['[', ']', ','], "");
                for word in stripped.split(' ') {
                    if let Some(id) = word.strip_prefix("UniMod:") {
                        let id = id
                            .parse::<i32>()
                            .map_err(|_| error(line, "invalid OBO UniMod record ID"))?;
                        term.record.record_id = u32::try_from(id).ok();
                    }
                }
            }
            "synonym" => {
                let synonym = text
                    .split('"')
                    .nth(1)
                    .filter(|_| text.matches('"').count() >= 2)
                    .ok_or_else(|| error(line, "OBO synonym requires quoted text"))?;
                synonym_count = checked(synonym_count, 1, options.max_aliases, "OBO synonyms")?;
                term.record.synonyms.insert(synonym.into());
            }
            "property_value" => {
                let value: String = value.chars().filter(|c| !c.is_whitespace()).collect();
                if value.contains("\"none\"") {
                    continue;
                }
                let parts: Vec<&str> = value.split('"').collect();
                if parts.len() != 3 {
                    return Err(error(line, "OBO property requires one quoted value"));
                }
                let property = parts[0];
                let argument = parts[1];
                let mass = || {
                    argument
                        .parse::<f64>()
                        .ok()
                        .filter(|v| v.is_finite())
                        .ok_or_else(|| error(line, "OBO mass must be finite"))
                };
                let formula =
                    || EmpiricalFormula::parse(argument).map_err(|e| error(line, e.to_string()));
                match property {
                    "DiffAvg:" => term.record.diff_average_mass = mass()?,
                    "DiffFormula:" => term.record.diff_formula = formula()?,
                    "DiffMono:" | "monoisotopicMass:" => term.record.diff_mono_mass = mass()?,
                    "Formula:" => {
                        term.record.absolute_formula = if argument.is_empty() {
                            None
                        } else {
                            Some(formula()?)
                        }
                    }
                    "MassAvg:" => term.record.average_mass = mass()?,
                    "MassMono:" => term.record.mono_mass = mass()?,
                    "Origin:" => term.origins = argument.into(),
                    "Source:" => term.record.classification = classification(argument).into(),
                    "TermSpec:" => {
                        term.record.term_specificity = match argument {
                            "C-term" => TermSpecificity::CTerm,
                            "N-term" => TermSpecificity::NTerm,
                            "none" => TermSpecificity::Anywhere,
                            // Source whitespace removal makes these otherwise valid
                            // ResidueModification spellings unreachable. Accept them.
                            "ProteinN-term" => TermSpecificity::ProteinNTerm,
                            "ProteinC-term" => TermSpecificity::ProteinCTerm,
                            _ => return Err(error(line, "invalid OBO terminal specificity")),
                        }
                    }
                    "reactionSites:" => {
                        if (argument == "2" && !options.cross_links_only)
                            || (argument == "1" && options.cross_links_only)
                        {
                            term.skip = true;
                        }
                    }
                    "specificities:" => {
                        term.origins = argument.replace(['(', ')'], "").replace('&', ",")
                    }
                    _ => {}
                }
            }
            _ => {}
        }
    }
    parsed.push_term(term, options)?;
    parsed.records.sort_by(|a, b| a.0.cmp(&b.0));
    Ok(parsed
        .records
        .into_iter()
        .map(|(_, record)| record)
        .collect())
}

fn read_line(
    reader: &mut impl BufRead,
    total: &mut usize,
    options: &OboReadOptions,
) -> Result<Option<Vec<u8>>> {
    let mut line = Vec::new();
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok((!line.is_empty()).then_some(line));
        }
        let length = buffer
            .iter()
            .position(|&byte| byte == b'\n')
            .map_or(buffer.len(), |index| index + 1);
        let done = buffer[length - 1] == b'\n';
        checked(line.len(), length, options.max_line_bytes, "OBO line bytes")?;
        *total = checked(*total, length, options.max_input_bytes, "OBO input bytes")?;
        line.extend_from_slice(&buffer[..length]);
        reader.consume(length);
        if done {
            return Ok(Some(line));
        }
    }
}

fn classification(value: &str) -> &'static str {
    match value.to_ascii_lowercase().as_str() {
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
