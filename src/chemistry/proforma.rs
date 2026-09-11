// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned ProForma annotation data, parsing, serialization and mass operations.
//!
//! This module implements the annotation AST, both text grammars, structured
//! errors, text writers and AASequence conversion from the pinned source.
//! Spectrum generation is not yet available. Resolution, conversion and mass/mz use
//! an explicitly supplied mutable registry. JSON
//! transport is available with the `proforma-json` feature.
//! Serialization preserves source omissions and does not validate ProForma
//! grammar. See `docs/PROFORMA_SUPPORT.md`, `docs/PROFORMA_PARSER_SUPPORT.md`
//! and `docs/PROFORMA_MASS_SUPPORT.md`.

use super::ResidueModification;
use crate::{Error, Result};
use std::sync::Arc;

/// Maximum serialized bytes per call, shared across every chain.
pub const MAX_PROFORMA_TEXT_BYTES: usize = 4 * 1024 * 1024;
/// Maximum consumed collection items per call.
pub const MAX_PROFORMA_NODES: usize = 1_000_000;
/// Maximum charged traversal, copying and numeric-formatting work per call.
pub const MAX_PROFORMA_WORK: usize = 50_000_000;

/// Source conversion policies. Both permissive variants currently have identical source behavior.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConversionPolicy {
    FailOnLoss,
    DropUnlocalised,
    BestEffort,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ConversionIssueType {
    UnresolvedMod,
    UnlocalisedMod,
    LabileMod,
    GlobalMod,
    AmbiguousMod,
    AmbiguousRegion,
    ModifiedRange,
    CrossLink,
    MultipleChains,
    AlternativeMods,
    UnsupportedFeature,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConversionIssue {
    pub issue_type: ConversionIssueType,
    pub description: String,
    /// None replaces the source's SIZE_MAX unknown-position sentinel.
    pub position: Option<usize>,
}

/// Both modes preserve vector order; source canonical writing does not sort.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum WriteMode {
    #[default]
    Lossless,
    Canonical,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CvDatabase {
    Unimod,
    Mod,
    Resid,
    Xlmod,
    Gno,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CvAccession {
    pub database: CvDatabase,
    pub accession: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NamedMod {
    pub cv_hint: Option<CvDatabase>,
    pub name: String,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MassDeltaSource {
    #[default]
    None,
    Obs,
    U,
    M,
    R,
    X,
    G,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct MassDelta {
    pub source: MassDeltaSource,
    pub mass: f64,
    /// Lossless mode copies a nonempty spelling without consulting mass.
    pub original_text: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FormulaTag {
    pub formula_string: String,
    pub charge: Option<i32>,
}

/// The source GlycanComposition::Monosaccharide variant, independent of the
/// immutable monosaccharide registry's record type.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GlycanComponent {
    Name(String),
    Formula(FormulaTag),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GlycanComposition {
    pub components: Vec<(GlycanComponent, i32)>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct InfoTag {
    pub text: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PositionConstraint {
    pub residues: Vec<char>,
    pub n_term: bool,
    pub c_term: bool,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ModificationTag {
    CvAccession(CvAccession),
    NamedMod(NamedMod),
    MassDelta(MassDelta),
    FormulaTag(FormulaTag),
    GlycanComposition(GlycanComposition),
    InfoTag(InfoTag),
    PositionConstraint(PositionConstraint),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LabelType {
    Crosslink,
    Branch,
    Ambiguous,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Label {
    /// Retained in the data model, but the source writer uses only identifier.
    pub label_type: LabelType,
    pub identifier: String,
    pub score: Option<f64>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Modification {
    pub alternatives: Vec<(ModificationTag, Option<Label>)>,
    /// Owned immutable chemistry, ignored by the text writer. Clones share it.
    pub resolved_mod: Option<Arc<ResidueModification>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SequenceElement {
    pub amino_acid: char,
    pub modifications: Vec<Modification>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct AmbiguousRegion {
    pub elements: Vec<SequenceElement>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ModifiedRange {
    pub elements: Vec<SequenceElement>,
    pub modifications: Vec<Modification>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SequenceSection {
    Element(SequenceElement),
    AmbiguousRegion(AmbiguousRegion),
    ModifiedRange(ModifiedRange),
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct UnlocalisedMod {
    pub modifications: Vec<Modification>,
    pub occurrence: Option<i32>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct LabileModification {
    pub modification: Modification,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct GlobalModification {
    pub modification: Modification,
    pub locations: Vec<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct IsotopeReplacement {
    pub isotope: String,
}

#[derive(Clone, Debug, PartialEq)]
pub enum GlobalModEntry {
    IsotopeReplacement(IsotopeReplacement),
    GlobalModification(GlobalModification),
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AdductIon {
    pub formula: String,
    pub charge: i32,
    pub occurrence: Option<i32>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChargeState {
    Simple(i32),
    Adducts(Vec<AdductIon>),
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Peptidoform {
    pub name: Option<String>,
    pub global_mods: Vec<GlobalModEntry>,
    pub unlocalised_mods: Vec<UnlocalisedMod>,
    pub labile_mods: Vec<LabileModification>,
    pub n_term_mods: Vec<Modification>,
    pub sequence: Vec<SequenceSection>,
    pub c_term_mods: Vec<Modification>,
    /// Source emits this only when the chain belongs to a chimeric ion.
    pub charge: Option<ChargeState>,
}

#[derive(Clone, Debug, Default, PartialEq)]
pub struct PeptidoformIon {
    /// Preserved in memory; source text serialization omits the ion name.
    pub name: Option<String>,
    pub chains: Vec<Peptidoform>,
    pub charge: Option<ChargeState>,
    pub is_chimeric: bool,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CrossLinkGroup {
    pub label: String,
    pub sites: Vec<(usize, usize)>,
}

impl Peptidoform {
    /// Serialize the complete source single-chain AST, omitting its charge.
    /// Errors publish no output and leave this value untouched.
    pub fn to_text(&self, mode: WriteMode) -> Result<String> {
        let mut writer = Writer::new(mode);
        writer.chain(self)?;
        Ok(writer.output)
    }
}

impl PeptidoformIon {
    /// Serialize all chains and source-supported charge placements atomically.
    /// The stored ion name is omitted, including in canonical mode.
    pub fn to_text(&self, mode: WriteMode) -> Result<String> {
        let mut writer = Writer::new(mode);
        writer.items(self.chains.len())?;
        for (i, chain) in self.chains.iter().enumerate() {
            if i != 0 {
                writer.text(if self.is_chimeric { "+" } else { "//" })?;
            }
            // Source serializes each chain in a fresh stream before insertion.
            writer.fixed_precision = None;
            writer.chain(chain)?;
            if self.is_chimeric {
                if let Some(charge) = &chain.charge {
                    writer.charge(charge)?;
                }
            }
        }
        if let Some(charge) = &self.charge {
            writer.charge(charge)?;
        }
        Ok(writer.output)
    }
}

fn invalid(message: &'static str) -> Error {
    Error::InvalidValue(message.into())
}

struct Writer {
    output: String,
    mode: WriteMode,
    fixed_precision: Option<usize>,
    work: usize,
    nodes: usize,
}

impl Writer {
    fn new(mode: WriteMode) -> Self {
        Self {
            output: String::new(),
            mode,
            fixed_precision: None,
            work: MAX_PROFORMA_WORK,
            nodes: MAX_PROFORMA_NODES,
        }
    }
    fn consume(&mut self, count: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(count)
            .ok_or_else(|| invalid("ProForma text work limit exceeded"))?;
        Ok(())
    }
    fn items(&mut self, count: usize) -> Result<()> {
        self.nodes = self
            .nodes
            .checked_sub(count)
            .ok_or_else(|| invalid("ProForma text node limit exceeded"))?;
        self.consume(count)
    }
    fn text(&mut self, value: &str) -> Result<()> {
        self.consume(value.len().saturating_add(1))?;
        let wanted = self
            .output
            .len()
            .checked_add(value.len())
            .filter(|n| *n <= MAX_PROFORMA_TEXT_BYTES)
            .ok_or_else(|| invalid("ProForma text output limit exceeded"))?;
        if wanted > self.output.capacity() {
            self.consume(self.output.len())?;
            // Geometric, capped capacity keeps copying and owned output bounded.
            let capacity = wanted
                .max(self.output.capacity().saturating_mul(2))
                .clamp(64, MAX_PROFORMA_TEXT_BYTES);
            self.output
                .try_reserve_exact(capacity - self.output.len())
                .map_err(|_| invalid("cannot allocate ProForma text"))?;
        }
        self.output.push_str(value);
        Ok(())
    }
    fn character(&mut self, value: char) -> Result<()> {
        if !value.is_ascii() {
            return Err(invalid("ProForma source character fields require ASCII"));
        }
        self.text(value.encode_utf8(&mut [0; 4]))
    }
    fn integer(&mut self, value: i64) -> Result<()> {
        self.consume(32)?;
        self.text(&value.to_string())
    }
    fn signed(&mut self, value: i32) -> Result<()> {
        if value >= 0 {
            self.text("+")?;
        }
        self.integer(i64::from(value))
    }
    fn float(&mut self, value: f64) -> Result<()> {
        if !value.is_finite() {
            return Err(invalid(
                "ProForma text requires finite consumed masses and scores",
            ));
        }
        // Fixed f64 output uses at most 315 bytes; charge all temporary scalar
        // formatting and copying before creating those bounded strings.
        self.consume(2048)?;
        let text = if let Some(precision) = self.fixed_precision {
            format!("{value:.precision$}")
        } else {
            crate::param::ParamValue::Float(value).to_stream_text()?
        };
        self.text(&text)
    }
    fn chain(&mut self, chain: &Peptidoform) -> Result<()> {
        self.consume(1)?;
        if let Some(name) = &chain.name {
            self.text("(>")?;
            self.text(name)?;
            self.text(")")?;
        }
        self.items(chain.global_mods.len())?;
        for entry in &chain.global_mods {
            self.text("<")?;
            match entry {
                GlobalModEntry::IsotopeReplacement(isotope) => self.text(&isotope.isotope)?,
                GlobalModEntry::GlobalModification(global) => {
                    self.modification(&global.modification, false)?;
                    self.items(global.locations.len())?;
                    for (i, location) in global.locations.iter().enumerate() {
                        self.text(if i == 0 { "@" } else { "," })?;
                        self.text(location)?;
                    }
                }
            }
            self.text(">")?;
        }
        self.items(chain.unlocalised_mods.len())?;
        for unlocalised in &chain.unlocalised_mods {
            self.modifications(&unlocalised.modifications)?;
            if let Some(occurrence) = unlocalised.occurrence {
                self.text("^")?;
                self.integer(i64::from(occurrence))?;
            }
            self.text("?")?;
        }
        self.items(chain.labile_mods.len())?;
        for labile in &chain.labile_mods {
            self.modification(&labile.modification, true)?;
        }
        self.modifications(&chain.n_term_mods)?;
        if !chain.n_term_mods.is_empty() {
            self.text("-")?;
        }
        self.items(chain.sequence.len())?;
        for section in &chain.sequence {
            match section {
                SequenceSection::Element(element) => self.element(element)?,
                SequenceSection::AmbiguousRegion(region) => {
                    self.text("(?")?;
                    self.elements(&region.elements)?;
                    self.text(")")?;
                }
                SequenceSection::ModifiedRange(range) => {
                    self.text("(")?;
                    self.elements(&range.elements)?;
                    self.text(")")?;
                    self.modifications(&range.modifications)?;
                }
            }
        }
        if !chain.c_term_mods.is_empty() {
            self.text("-")?;
        }
        self.modifications(&chain.c_term_mods)
    }
    fn elements(&mut self, elements: &[SequenceElement]) -> Result<()> {
        self.items(elements.len())?;
        for element in elements {
            self.element(element)?;
        }
        Ok(())
    }
    fn element(&mut self, element: &SequenceElement) -> Result<()> {
        self.character(element.amino_acid)?;
        self.modifications(&element.modifications)
    }
    fn modifications(&mut self, modifications: &[Modification]) -> Result<()> {
        self.items(modifications.len())?;
        for modification in modifications {
            self.modification(modification, false)?;
        }
        Ok(())
    }
    fn modification(&mut self, modification: &Modification, labile: bool) -> Result<()> {
        self.items(modification.alternatives.len())?;
        self.text(if labile { "{" } else { "[" })?;
        for (i, (tag, label)) in modification.alternatives.iter().enumerate() {
            if i != 0 {
                self.text("|")?;
            }
            let label_only = !labile
                && label.is_some()
                && matches!(tag, ModificationTag::InfoTag(info) if info.text.is_empty());
            if !label_only {
                self.tag(tag)?;
            }
            if let Some(label) = label {
                self.text("#")?;
                self.text(&label.identifier)?;
                if let Some(score) = label.score {
                    self.text("(")?;
                    if self.mode == WriteMode::Canonical {
                        self.fixed_precision = Some(2);
                    }
                    self.float(score)?;
                    self.text(")")?;
                }
            }
        }
        self.text(if labile { "}" } else { "]" })
    }
    fn tag(&mut self, tag: &ModificationTag) -> Result<()> {
        match tag {
            ModificationTag::CvAccession(cv) => {
                self.text(match cv.database {
                    CvDatabase::Unimod => "UNIMOD",
                    CvDatabase::Mod => "MOD",
                    CvDatabase::Resid => "RESID",
                    CvDatabase::Xlmod => "XLMOD",
                    CvDatabase::Gno => "GNO",
                })?;
                self.text(":")?;
                self.text(&cv.accession)
            }
            ModificationTag::NamedMod(named) => {
                if let Some(hint) = named.cv_hint {
                    self.text(match hint {
                        CvDatabase::Unimod => "U:",
                        CvDatabase::Mod => "M:",
                        CvDatabase::Resid => "R:",
                        CvDatabase::Xlmod => "X:",
                        CvDatabase::Gno => "G:",
                    })?;
                }
                self.text(&named.name)
            }
            ModificationTag::MassDelta(delta) => {
                self.text(match delta.source {
                    MassDeltaSource::None => "",
                    MassDeltaSource::Obs => "Obs:",
                    MassDeltaSource::U => "U:",
                    MassDeltaSource::M => "M:",
                    MassDeltaSource::R => "R:",
                    MassDeltaSource::X => "X:",
                    MassDeltaSource::G => "G:",
                })?;
                if self.mode == WriteMode::Lossless && !delta.original_text.is_empty() {
                    self.text(&delta.original_text)
                } else {
                    if delta.mass >= 0.0 {
                        self.text("+")?;
                    }
                    self.fixed_precision = Some(4);
                    self.float(delta.mass)
                }
            }
            ModificationTag::FormulaTag(formula) => self.formula(formula),
            ModificationTag::GlycanComposition(glycan) => {
                self.items(glycan.components.len())?;
                self.text("Glycan:")?;
                for (component, count) in &glycan.components {
                    match component {
                        GlycanComponent::Name(name) => self.text(name)?,
                        GlycanComponent::Formula(formula) => self.formula(formula)?,
                    }
                    self.integer(i64::from(*count))?;
                }
                Ok(())
            }
            ModificationTag::InfoTag(info) => {
                self.text("INFO:")?;
                self.text(&info.text)
            }
            ModificationTag::PositionConstraint(position) => {
                self.items(position.residues.len())?;
                self.text("Position:")?;
                if position.n_term {
                    self.text("N-term")?;
                }
                if position.c_term {
                    if position.n_term {
                        self.text(",")?;
                    }
                    self.text("C-term")?;
                }
                if !position.residues.is_empty() && (position.n_term || position.c_term) {
                    self.text(",")?;
                }
                for residue in &position.residues {
                    self.character(*residue)?;
                }
                Ok(())
            }
        }
    }
    fn formula(&mut self, formula: &FormulaTag) -> Result<()> {
        self.text("Formula:")?;
        self.text(&formula.formula_string)?;
        if let Some(charge) = formula.charge {
            self.text(":z")?;
            self.signed(charge)?;
        }
        Ok(())
    }
    fn charge(&mut self, charge: &ChargeState) -> Result<()> {
        self.text("/")?;
        match charge {
            ChargeState::Simple(value) => self.integer(i64::from(*value)),
            ChargeState::Adducts(adducts) => {
                self.items(adducts.len())?;
                self.text("[")?;
                let mut total = 0i32;
                for (i, adduct) in adducts.iter().enumerate() {
                    if i != 0 {
                        self.text(",")?;
                    }
                    self.text(&adduct.formula)?;
                    self.text(":z")?;
                    self.signed(adduct.charge)?;
                    if let Some(occurrence) = adduct.occurrence {
                        self.text("^")?;
                        self.integer(i64::from(occurrence))?;
                    }
                    let contribution = adduct
                        .charge
                        .checked_mul(adduct.occurrence.unwrap_or(1))
                        .ok_or_else(|| invalid("ProForma adduct charge multiplication overflow"))?;
                    total = total
                        .checked_add(contribution)
                        .ok_or_else(|| invalid("ProForma adduct charge sum overflow"))?;
                }
                self.text("]")?;
                self.integer(i64::from(total.unsigned_abs()))?;
                self.text(if total >= 0 { "+" } else { "-" })
            }
        }
    }
}

#[path = "proforma_parser.rs"]
mod parser;
pub use parser::{
    ErrorCode, MAX_PROFORMA_DIAGNOSTIC_BYTES, MAX_PROFORMA_PARSE_BYTES, MAX_PROFORMA_PARSE_DEPTH,
    MAX_PROFORMA_PARSE_INPUT_BYTES, MAX_PROFORMA_PARSE_NODES, MAX_PROFORMA_PARSE_WORK, ParseError,
    ParseFailure,
};

#[path = "proforma_resolution.rs"]
mod resolution;
pub use resolution::{
    MAX_PROFORMA_RESOLUTION_BYTES, MAX_PROFORMA_RESOLUTION_ITEMS,
    MAX_PROFORMA_RESOLUTION_TEXT_BYTES, MAX_PROFORMA_RESOLUTION_WORK, ResolutionWarning,
};

#[cfg(feature = "proforma-json")]
#[path = "proforma_json.rs"]
mod json;
#[cfg(feature = "proforma-json")]
pub use json::{
    MAX_PROFORMA_JSON_BYTES, MAX_PROFORMA_JSON_DEPTH, MAX_PROFORMA_JSON_ITEMS,
    MAX_PROFORMA_JSON_TEXT_BYTES, MAX_PROFORMA_JSON_WORK,
};

#[path = "proforma_mass.rs"]
mod mass;
pub use mass::{
    MAX_PROFORMA_MASS_BYTES, MAX_PROFORMA_MASS_ITEMS, MAX_PROFORMA_MASS_TEXT_BYTES,
    MAX_PROFORMA_MASS_WORK, MassAttempt, MassEvaluation,
};

#[path = "proforma_conversion.rs"]
mod conversion;
pub use conversion::{
    ConversionEvaluation, ConversionWarning, MAX_PROFORMA_CONVERSION_BYTES,
    MAX_PROFORMA_CONVERSION_ITEMS, MAX_PROFORMA_CONVERSION_TEXT_BYTES,
    MAX_PROFORMA_CONVERSION_WORK,
};
