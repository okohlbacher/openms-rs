// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{
    EmpiricalFormula, FragmentIon, IonSeries, ModificationsDB, PROTON_MASS_U, ResidueModification,
    TermSpecificity, composition_formula, ion_mz, residue_composition,
};
use crate::{Error, Result};
use std::{cmp::Ordering, fmt, ops::Range, str::FromStr, sync::Arc};

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn unknown(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn valid_mass(mass: f64) -> Result<f64> {
    if mass.is_finite() && mass >= 0.0 {
        Ok(mass)
    } else {
        Err(invalid(
            "calculated sequence/residue mass must be finite and nonnegative",
        ))
    }
}
fn water() -> EmpiricalFormula {
    composition_formula([0, 2, 0, 1, 0, 0])
}
fn terminal_base(term: TermSpecificity) -> f64 {
    if matches!(term, TermSpecificity::NTerm | TermSpecificity::ProteinNTerm) {
        composition_formula([0, 1, 0, 0, 0, 0]).mono_mass()
    } else {
        composition_formula([0, 1, 0, 1, 0, 0]).mono_mass()
    }
}

/// All finite source peptide fragment types. Formula and mass queries retain
/// the source's internal-residue fallback for radical, precursor, neutral-loss
/// and unassigned variants; they do not generate a spectrum.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PeptideFragmentType {
    #[default]
    Full,
    Internal,
    NTerminal,
    CTerminal,
    AIon,
    BIon,
    CIon,
    XIon,
    YIon,
    ZIon,
    Zp1Ion,
    Zp2Ion,
    Precursor,
    BIonMinusH2O,
    YIonMinusH2O,
    BIonMinusNH3,
    YIonMinusNH3,
    NonIdentified,
    Unannotated,
}
impl PeptideFragmentType {
    fn keeps_n_terminal(self) -> bool {
        matches!(
            self,
            Self::Full | Self::NTerminal | Self::AIon | Self::BIon | Self::CIon
        )
    }
    fn keeps_c_terminal(self) -> bool {
        matches!(
            self,
            Self::Full | Self::CTerminal | Self::XIon | Self::YIon | Self::ZIon
        )
    }
    fn correction(self) -> Option<[i32; 6]> {
        match self {
            Self::Full | Self::YIon => Some([0, 2, 0, 1, 0, 0]),
            Self::NTerminal => Some([0, 1, 0, 0, 0, 0]),
            Self::CTerminal => Some([0, 1, 0, 1, 0, 0]),
            Self::AIon => Some([-1, 0, 0, -1, 0, 0]),
            Self::CIon => Some([0, 3, 1, 0, 0, 0]),
            Self::XIon => Some([1, 0, 0, 2, 0, 0]),
            Self::ZIon => Some([0, -1, -1, 1, 0, 0]),
            _ => None,
        }
    }
}

/// Owned anonymous monoisotopic annotation. The original decimal spelling is
/// preserved; no atom composition or average mass is inferred from one mass.
/// Constructed by AASequence parsing or its checked mass-tag setters.
#[derive(Clone, Debug, PartialEq)]
pub struct MassTag {
    // Fragment slices own references to immutable spellings, like known
    // modifications, rather than copying the same decimal token per fragment.
    input: Arc<str>,
    full_id: Arc<str>,
    mass: f64,
    is_delta: bool,
    delta_mono_mass: Option<f64>,
    residue_mono_mass: Option<f64>,
    origin: Option<char>,
    term: TermSpecificity,
}
// Private fields are constructed only from finite checked values.
impl Eq for MassTag {}
// Chemical values are finite. partial_cmp preserves Eq's signed-zero semantics.
impl Ord for MassTag {
    fn cmp(&self, other: &Self) -> Ordering {
        (
            &self.input,
            &self.full_id,
            self.is_delta,
            self.origin,
            self.term,
        )
            .cmp(&(
                &other.input,
                &other.full_id,
                other.is_delta,
                other.origin,
                other.term,
            ))
            .then_with(|| {
                self.mass
                    .partial_cmp(&other.mass)
                    .expect("validated finite mass")
            })
            .then_with(|| {
                self.delta_mono_mass
                    .partial_cmp(&other.delta_mono_mass)
                    .expect("validated finite mass")
            })
            .then_with(|| {
                self.residue_mono_mass
                    .partial_cmp(&other.residue_mono_mass)
                    .expect("validated finite mass")
            })
    }
}
impl PartialOrd for MassTag {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl MassTag {
    pub fn full_id(&self) -> &str {
        &self.full_id
    }
    pub fn origin(&self) -> Option<char> {
        self.origin
    }
    pub fn term_specificity(&self) -> TermSpecificity {
        self.term
    }
    pub fn input(&self) -> &str {
        &self.input
    }
    /// Number written inside brackets; a sign denotes a mass difference.
    pub fn mass(&self) -> f64 {
        self.mass
    }
    pub fn is_delta(&self) -> bool {
        self.is_delta
    }
    /// None for an absolute tag on B/Z/X, whose unmodified mass is unknown.
    pub fn delta_mono_mass(&self) -> Option<f64> {
        self.delta_mono_mass
    }
    /// Internal residue mass, including its tag; None for terminal annotations.
    pub fn residue_mono_mass(&self) -> Option<f64> {
        self.residue_mono_mass
    }
}

/// An immutable registry modification or an owned anonymous mass annotation.
/// Known chemistry is shared and owned; parsing never mutates the database.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum SequenceModification {
    Known(Arc<ResidueModification>),
    MassTag(MassTag),
}
impl SequenceModification {
    pub fn known(&self) -> Option<&ResidueModification> {
        match self {
            Self::Known(m) => Some(m.as_ref()),
            Self::MassTag(_) => None,
        }
    }
    pub fn mass_tag(&self) -> Option<&MassTag> {
        match self {
            Self::Known(_) => None,
            Self::MassTag(tag) => Some(tag),
        }
    }
    pub fn name(&self) -> &str {
        match self {
            Self::Known(m) => m.name(),
            Self::MassTag(tag) => &tag.input,
        }
    }
    pub fn full_id(&self) -> &str {
        match self {
            Self::Known(m) => m.full_id(),
            Self::MassTag(tag) => &tag.full_id,
        }
    }
    pub fn record_id(&self) -> Option<u32> {
        self.known().and_then(ResidueModification::record_id)
    }
    pub fn origin(&self) -> Option<char> {
        match self {
            Self::Known(m) => m.origin(),
            Self::MassTag(tag) => tag.origin,
        }
    }
    pub fn term_specificity(&self) -> TermSpecificity {
        match self {
            Self::Known(m) => m.term_specificity(),
            Self::MassTag(tag) => tag.term,
        }
    }
    pub fn diff_mono_mass(&self) -> Result<f64> {
        match self {
            Self::Known(m) => Ok(m.diff_mono_mass()),
            Self::MassTag(tag) => tag.delta_mono_mass.ok_or_else(|| {
                unknown("absolute tag has no known unmodified residue mass difference")
            }),
        }
    }
    pub fn diff_formula(&self) -> Result<&EmpiricalFormula> {
        self.known()
            .filter(|m| !missing_modification_formula(m))
            .map(ResidueModification::diff_formula)
            .ok_or_else(|| unknown("mass annotation has no known atom formula"))
    }
    pub fn diff_average_mass(&self) -> Result<f64> {
        self.known()
            .map(ResidueModification::diff_average_mass)
            .ok_or_else(|| unknown("anonymous monoisotopic annotation has no known average mass"))
    }
}

/// Sequence with residue/terminal annotations and independently known chemistry.
/// B/Z/X are retained for identification and digestion; unresolved chemistry
/// returns Error::Unsupported. Numeric brackets resolve known modifications at
/// the source's written precision, or retain an owned mass tag. See
/// `docs/SEQUENCE_SUPPORT.md` for source conventions and native corrections.
///
/// ```
/// use openms::chemistry::AASequence;
/// let peptide = AASequence::parse("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K")?;
/// assert_eq!(peptide.as_str(), "ACMK");
/// assert_eq!(AASequence::parse(&peptide.to_string())?, peptide);
/// assert_eq!(peptide.fragment_ions(2)?.len(), 12);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Clone, Debug, PartialEq)]
pub struct AASequence {
    sequence: String,
    formula: Option<EmpiricalFormula>,
    mono_mass: Option<f64>,
    residue_modifications: Vec<Option<SequenceModification>>,
    n_terminal: Option<SequenceModification>,
    c_terminal: Option<SequenceModification>,
}
impl Eq for AASequence {}
// Deterministic native value order. Chemistry is part of identity even when two
// caller registries use the same displayed identifiers for different records.
impl Ord for AASequence {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sequence
            .cmp(&other.sequence)
            .then_with(|| self.n_terminal.cmp(&other.n_terminal))
            .then_with(|| self.residue_modifications.cmp(&other.residue_modifications))
            .then_with(|| self.c_terminal.cmp(&other.c_terminal))
            .then_with(|| {
                self.formula
                    .as_ref()
                    .map(|formula| (&formula.atoms, formula.charge))
                    .cmp(
                        &other
                            .formula
                            .as_ref()
                            .map(|formula| (&formula.atoms, formula.charge)),
                    )
            })
            .then_with(|| {
                self.mono_mass
                    .partial_cmp(&other.mono_mass)
                    .expect("validated finite sequence mass")
            })
    }
}
impl PartialOrd for AASequence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Default for AASequence {
    fn default() -> Self {
        Self {
            sequence: String::new(),
            formula: Some(EmpiricalFormula::default()),
            mono_mass: Some(0.0),
            residue_modifications: Vec::new(),
            n_terminal: None,
            c_terminal: None,
        }
    }
}
impl AASequence {
    pub fn parse(input: &str) -> Result<Self> {
        input.parse()
    }
    /// Unmodified uppercase residue string; Display includes annotations.
    pub fn as_str(&self) -> &str {
        &self.sequence
    }
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
    /// Neutral full formula, if every residue and annotation has known atoms.
    pub fn formula(&self) -> Result<EmpiricalFormula> {
        self.formula.clone().ok_or_else(|| {
            unknown("sequence contains unresolved residues or annotations without an atom formula")
        })
    }
    /// Formula of this entire sequence treated as the requested fragment type.
    /// Slice first to select residues. Charge is formula metadata, not added H
    /// atoms. Empty sequences ignore both charge and annotations. Only retained
    /// termini require known composition; signed fragment atom counts are valid.
    /// Each call permits 50 million work units and 256 MiB cumulative formula
    /// allocation accounting, including temporary maps.
    pub fn formula_for(
        &self,
        fragment: PeptideFragmentType,
        charge: i32,
    ) -> Result<EmpiricalFormula> {
        self.formula_for_with_budget(fragment, charge, &mut 50_000_000, &mut (256 * 1024 * 1024))
    }

    /// Shared graph-query budget. Charges survive errors as well as success.
    pub(crate) fn formula_for_with_budget(
        &self,
        fragment: PeptideFragmentType,
        charge: i32,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<EmpiricalFormula> {
        if self.is_empty() {
            return Ok(EmpiricalFormula::default());
        }
        *remaining_work = remaining_work
            .checked_sub(self.len())
            .ok_or_else(|| invalid("peptide formula work limit exceeded"))?;
        let mut result = EmpiricalFormula::default().with_charge(charge);
        // AASequence.cpp adds relevant N/C deltas before internal residues.
        for modification in [
            self.n_terminal
                .as_ref()
                .filter(|_| fragment.keeps_n_terminal()),
            self.c_terminal
                .as_ref()
                .filter(|_| fragment.keeps_c_terminal()),
        ]
        .into_iter()
        .flatten()
        {
            add_formula(
                &mut result,
                modification.diff_formula()?,
                remaining_work,
                remaining_bytes,
            )?;
        }
        for index in 0..self.len() {
            // Precharge internal_formula's base, delta/replacement and water
            // maps before it clones any caller-provided chemical formula.
            reserve_formula(6, remaining_work, remaining_bytes)?;
            if let Some(SequenceModification::Known(modification)) =
                &self.residue_modifications[index]
            {
                if !modification.diff_formula().is_empty() {
                    reserve_formula(
                        6usize.saturating_add(modification.diff_formula().atoms.len()),
                        remaining_work,
                        remaining_bytes,
                    )?;
                } else if missing_modification_formula(modification) {
                    if let Some(absolute) = modification.absolute_formula() {
                        reserve_formula(2, remaining_work, remaining_bytes)?;
                        reserve_formula(
                            absolute.atoms.len().saturating_add(2),
                            remaining_work,
                            remaining_bytes,
                        )?;
                    }
                }
            }
            let internal = self.internal_formula(index)?.ok_or_else(|| {
                unknown(
                    "fragment contains an unresolved residue or annotation without an atom formula",
                )
            })?;
            add_formula(&mut result, &internal, remaining_work, remaining_bytes)?;
        }
        let Some(correction) = fragment.correction() else {
            return Ok(result);
        };
        reserve_formula(3, remaining_work, remaining_bytes)?;
        add_formula(
            &mut result,
            &composition_formula(correction),
            remaining_work,
            remaining_bytes,
        )?;
        Ok(result)
    }
    /// Full neutral monoisotopic mass. Known residue modifications prefer formula
    /// masses, with source declared-mass fallback when their formula is absent;
    /// known terminal modifications retain declared source delta masses.
    pub fn mono_mass(&self) -> Result<f64> {
        self.mono_mass
            .ok_or_else(|| unknown("sequence contains a residue without a known monoisotopic mass"))
    }
    /// Monoisotopic ion mass for all supplied residues and a source fragment
    /// type. Includes `charge * PROTON_MASS_U`, without dividing by charge.
    /// Source scalar order is charge, N/C deltas, residues, then correction;
    /// rounding can differ from the cached neutral `mono_mass()` calculation.
    /// Mass-only annotations are supported, and finite signed outputs are kept.
    pub fn mono_mass_for(&self, fragment: PeptideFragmentType, charge: i32) -> Result<f64> {
        self.mono_mass_for_with_budget(fragment, charge, &mut 50_000_000, &mut (256 * 1024 * 1024))
    }
    fn mono_mass_for_with_budget(
        &self,
        fragment: PeptideFragmentType,
        charge: i32,
        work: &mut usize,
        bytes: &mut usize,
    ) -> Result<f64> {
        if self.is_empty() {
            return Ok(0.0);
        }
        consume_mass_work(self.len(), work)?;
        let mut mass = f64::from(charge) * PROTON_MASS_U;
        for modification in [
            self.n_terminal
                .as_ref()
                .filter(|_| fragment.keeps_n_terminal()),
            self.c_terminal
                .as_ref()
                .filter(|_| fragment.keeps_c_terminal()),
        ]
        .into_iter()
        .flatten()
        {
            mass = finite_fragment_mass(mass + modification.diff_mono_mass()?)?;
        }
        for index in 0..self.len() {
            // internal_mass builds the ordinary base formula even for a mass
            // tag. Formula-free records can additionally build full-residue
            // water/base maps; count them before calling the shared chemistry.
            reserve_formula(6, work, bytes)?;
            if let Some(SequenceModification::Known(modification)) =
                &self.residue_modifications[index]
            {
                if missing_modification_formula(modification) {
                    reserve_formula(2, work, bytes)?;
                    if let Some(absolute) = modification.absolute_formula() {
                        consume_mass_work(absolute.atoms.len(), work)?;
                    } else if modification.mono_mass() == 0.0 {
                        reserve_formula(6, work, bytes)?;
                        reserve_formula(8, work, bytes)?;
                    }
                } else {
                    consume_mass_work(modification.diff_formula().atoms.len(), work)?;
                }
            }
            let internal = self.internal_mass(index)?.ok_or_else(|| {
                unknown("fragment contains a residue without a known monoisotopic mass")
            })?;
            mass = finite_fragment_mass(mass + internal)?;
        }
        if let Some(correction) = fragment.correction() {
            reserve_formula(3, work, bytes)?;
            mass = finite_fragment_mass(mass + composition_formula(correction).mono_mass())?;
        }
        Ok(mass)
    }
    /// Abundance-weighted average from the complete known formula.
    pub fn average_mass(&self) -> Result<f64> {
        self.formula
            .as_ref()
            .map(EmpiricalFormula::average_mass)
            .ok_or_else(|| unknown("sequence has no complete formula for an average mass"))
    }
    /// Abundance-weighted ion mass from the requested fragment formula. No
    /// average mass is fabricated for anonymous or formula-free annotations.
    /// Uses the same terminal selection, charge and limits as `formula_for`.
    pub fn average_mass_for(&self, fragment: PeptideFragmentType, charge: i32) -> Result<f64> {
        let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
        let formula = self.formula_for_with_budget(fragment, charge, &mut work, &mut bytes)?;
        consume_mass_work(formula.atoms.len(), &mut work)?;
        finite_fragment_mass(formula.average_mass())
    }
    pub fn mz(&self, charge: i32) -> Result<f64> {
        if self.is_empty() {
            return Err(invalid("an empty sequence has no ion m/z"));
        }
        ion_mz(
            self.mono_mass()? + f64::from(charge) * PROTON_MASS_U,
            charge,
        )
    }
    pub fn is_modified(&self) -> bool {
        self.n_terminal.is_some()
            || self.c_terminal.is_some()
            || self.residue_modifications.iter().any(Option::is_some)
    }
    pub fn n_terminal_modification(&self) -> Option<&SequenceModification> {
        self.n_terminal.as_ref()
    }
    pub fn c_terminal_modification(&self) -> Option<&SequenceModification> {
        self.c_terminal.as_ref()
    }
    pub fn residue_modification(&self, index: usize) -> Result<Option<&SequenceModification>> {
        self.residue_modifications
            .get(index)
            .map(Option::as_ref)
            .ok_or_else(|| invalid("residue index out of bounds"))
    }
    fn residue(&self, index: usize) -> Result<char> {
        self.sequence
            .as_bytes()
            .get(index)
            .copied()
            .map(char::from)
            .ok_or_else(|| invalid("residue index out of bounds"))
    }
    /// Set/replace a named residue modification; empty name removes it.
    /// Any error leaves the sequence unchanged.
    pub fn set_modification(&mut self, index: usize, name: &str) -> Result<()> {
        self.set_modification_with_registry(index, name, ModificationsDB::global())
    }
    /// Resolve against a caller-owned registry; the peptide retains shared chemistry.
    pub fn set_modification_with_registry(
        &mut self,
        index: usize,
        name: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        let residue = self.residue(index)?;
        let modification = if name.is_empty() {
            None
        } else {
            Some(SequenceModification::Known(db.get_modification_handle(
                name,
                Some(residue),
                Some(TermSpecificity::Anywhere),
            )?))
        };
        self.replace(Attachment::Residue(index), modification)
    }
    /// Set a numeric residue tag, using the same precision/lookup as brackets.
    /// `text` is the bracket content. This replaces the existing annotation.
    pub fn set_mass_tag(&mut self, index: usize, text: &str) -> Result<()> {
        self.set_mass_tag_with_registry(index, text, ModificationsDB::global())
    }
    /// Resolve a numeric tag against a supplied registry, retaining unknown tags.
    pub fn set_mass_tag_with_registry(
        &mut self,
        index: usize,
        text: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        let residue = self.residue(index)?;
        let modification = resolve_mass(text, residue, TermSpecificity::Anywhere, db)?;
        self.replace(Attachment::Residue(index), Some(modification))
    }
    pub fn set_n_terminal_modification(&mut self, name: &str) -> Result<()> {
        self.set_n_terminal_modification_with_registry(name, ModificationsDB::global())
    }
    pub fn set_c_terminal_modification(&mut self, name: &str) -> Result<()> {
        self.set_c_terminal_modification_with_registry(name, ModificationsDB::global())
    }
    pub fn set_n_terminal_mass_tag(&mut self, text: &str) -> Result<()> {
        self.set_n_terminal_mass_tag_with_registry(text, ModificationsDB::global())
    }
    pub fn set_c_terminal_mass_tag(&mut self, text: &str) -> Result<()> {
        self.set_c_terminal_mass_tag_with_registry(text, ModificationsDB::global())
    }
    pub fn set_n_terminal_modification_with_registry(
        &mut self,
        name: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        self.set_terminal(name, true, false, db)
    }
    pub fn set_c_terminal_modification_with_registry(
        &mut self,
        name: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        self.set_terminal(name, false, false, db)
    }
    pub fn set_n_terminal_mass_tag_with_registry(
        &mut self,
        text: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        self.set_terminal(text, true, true, db)
    }
    pub fn set_c_terminal_mass_tag_with_registry(
        &mut self,
        text: &str,
        db: &ModificationsDB,
    ) -> Result<()> {
        self.set_terminal(text, false, true, db)
    }
    fn set_terminal(
        &mut self,
        text: &str,
        n_terminal: bool,
        mass: bool,
        db: &ModificationsDB,
    ) -> Result<()> {
        if self.is_empty() {
            return Err(invalid("empty sequence has no terminus"));
        }
        let residue = self.residue(if n_terminal { 0 } else { self.len() - 1 })?;
        let modification = if !mass && text.is_empty() {
            None
        } else if mass {
            Some(resolve_mass(
                text,
                residue,
                if n_terminal {
                    TermSpecificity::NTerm
                } else {
                    TermSpecificity::CTerm
                },
                db,
            )?)
        } else {
            Some(SequenceModification::Known(resolve_terminal(
                text, residue, n_terminal, db,
            )?))
        };
        self.replace(
            if n_terminal {
                Attachment::NTerm
            } else {
                Attachment::CTerm
            },
            modification,
        )
    }
    fn replace(
        &mut self,
        attachment: Attachment,
        modification: Option<SequenceModification>,
    ) -> Result<()> {
        let mut next = self.clone();
        *next.slot(attachment) = modification;
        next.rebuild_chemistry()?;
        *self = next;
        Ok(())
    }
    // The source generator installs resolved residue pointers directly, including
    // terminal-specific records on residue slots in its maximum-one fast path.
    // Keep this escape hatch internal; ordinary public setters still check origin.
    pub(crate) fn with_resolved_modifications(
        &self,
        residues: &[(usize, Arc<ResidueModification>)],
        n_terminal: Option<Arc<ResidueModification>>,
        c_terminal: Option<Arc<ResidueModification>>,
    ) -> Result<Self> {
        if residues.iter().any(|(index, _)| *index >= self.len()) {
            return Err(invalid("resolved modification index out of bounds"));
        }
        let mut next = self.clone();
        for (index, modification) in residues {
            next.residue_modifications[*index] =
                Some(SequenceModification::Known(Arc::clone(modification)));
        }
        if let Some(modification) = n_terminal {
            next.n_terminal = Some(SequenceModification::Known(modification));
        }
        if let Some(modification) = c_terminal {
            next.c_terminal = Some(SequenceModification::Known(modification));
        }
        if next.is_empty() {
            // Source empty sequences retain typed terminal pointers but report
            // zero chemistry, independently of those annotations.
            next.formula = Some(EmpiricalFormula::default());
            next.mono_mass = Some(0.0);
        } else {
            next.rebuild_chemistry()?;
        }
        Ok(next)
    }

    // Logical owned-payload allowance for generation, without constructing a
    // display string or cloning unknown tags. Formula entries are deliberately
    // overcounted: six standard-residue element types plus every known delta.
    pub(crate) fn generation_payload_bytes(&self) -> Result<usize> {
        let overflow = || invalid("sequence payload size overflows");
        let mut bytes = std::mem::size_of::<Self>()
            .checked_add(self.sequence.len())
            .and_then(|n| {
                self.residue_modifications
                    .len()
                    .checked_mul(std::mem::size_of::<Option<SequenceModification>>())
                    .and_then(|slots| n.checked_add(slots))
            })
            .ok_or_else(overflow)?;
        let mut formula_entries = 6usize;
        for modification in self
            .residue_modifications
            .iter()
            .chain([&self.n_terminal, &self.c_terminal])
            .flatten()
        {
            match modification {
                SequenceModification::Known(m) => {
                    formula_entries = formula_entries
                        .checked_add(m.diff_formula().atoms.len())
                        .and_then(|n| {
                            n.checked_add(m.absolute_formula().map_or(0, |f| f.atoms.len()))
                        })
                        .ok_or_else(overflow)?;
                }
                SequenceModification::MassTag(tag) => {
                    bytes = bytes
                        .checked_add(tag.input.len())
                        .and_then(|n| n.checked_add(tag.full_id.len()))
                        .ok_or_else(overflow)?;
                }
            }
        }
        formula_entries
            .checked_mul(GENERATION_FORMULA_ENTRY_BYTES)
            .and_then(|n| bytes.checked_add(n))
            .ok_or_else(overflow)
    }
    fn slot(&mut self, attachment: Attachment) -> &mut Option<SequenceModification> {
        match attachment {
            Attachment::NTerm => &mut self.n_terminal,
            Attachment::CTerm => &mut self.c_terminal,
            Attachment::Residue(index) => &mut self.residue_modifications[index],
        }
    }
    /// Slice by residue indices. End modifications survive only when their
    /// original terminus survives. Known chemistry is recomputed for the slice.
    pub fn subsequence(&self, range: Range<usize>) -> Result<Self> {
        let sequence = self
            .sequence
            .get(range.clone())
            .ok_or_else(|| invalid("peptide subsequence range is out of bounds"))?;
        if sequence.is_empty() {
            return Ok(Self::default());
        }
        let mut result = Self {
            sequence: sequence.into(),
            residue_modifications: self.residue_modifications[range.clone()].to_vec(),
            n_terminal: if range.start == 0 {
                self.n_terminal.clone()
            } else {
                None
            },
            c_terminal: if range.end == self.len() {
                self.c_terminal.clone()
            } else {
                None
            },
            ..Self::default()
        };
        result.rebuild_chemistry()?;
        Ok(result)
    }
    pub fn prefix(&self, length: usize) -> Result<Self> {
        self.subsequence(0..length)
    }
    pub fn suffix(&self, length: usize) -> Result<Self> {
        self.subsequence(
            self.len()
                .checked_sub(length)
                .ok_or_else(|| invalid("suffix length out of bounds"))?..self.len(),
        )
    }
    /// Render UniMod accessions where available, otherwise absolute mass brackets.
    /// Anonymous tags keep their original spelling. Non-UniMod names/chemistry are
    /// intentionally lost in this mass representation; unknown masses are errors.
    pub fn to_unimod_string(&self) -> Result<String> {
        if self.is_empty() {
            // Source exports an empty string even when its generator attached
            // typed terminal records to an empty sequence.
            return Ok(String::new());
        }
        let annotation = |m: &SequenceModification,
                          internal: Option<f64>,
                          term: Option<TermSpecificity>|
         -> Result<String> {
            if let SequenceModification::MassTag(tag) = m {
                return Ok(format!("[{}]", tag.input));
            }
            if let Some(accession) = m.known().and_then(ResidueModification::unimod_accession) {
                return Ok(format!("({accession})"));
            }
            let mass = if let Some(term) = term {
                valid_mass(terminal_base(term) + m.diff_mono_mass()?)?
            } else {
                internal.ok_or_else(|| {
                    unknown("non-UniMod residue annotation has no known absolute mass")
                })?
            };
            Ok(format!("[{mass}]"))
        };
        let mut result = String::new();
        if let Some(m) = &self.n_terminal {
            result.push('.');
            result.push_str(&annotation(m, None, Some(TermSpecificity::NTerm))?);
        }
        for (index, (residue, modification)) in self
            .sequence
            .chars()
            .zip(&self.residue_modifications)
            .enumerate()
        {
            result.push(residue);
            if let Some(m) = modification {
                result.push_str(&annotation(m, self.internal_mass(index)?, None)?);
            }
        }
        if let Some(m) = &self.c_terminal {
            result.push('.');
            result.push_str(&annotation(m, None, Some(TermSpecificity::CTerm))?);
        }
        Ok(result)
    }
    /// Render each known vocabulary accession, retaining custom names when no
    /// accession exists. Exact reconstruction requires the same registry.
    pub fn to_accession_string(&self) -> String {
        self.annotated_string(true)
    }
    fn annotated_string(&self, accessions: bool) -> String {
        let annotation = |m: &SequenceModification| match m {
            SequenceModification::MassTag(tag) => format!("[{}]", tag.input),
            SequenceModification::Known(m) => {
                let name = if accessions && !m.accession().is_empty() {
                    m.accession()
                } else if m.name().is_empty()
                    || matches!(
                        m.term_specificity(),
                        TermSpecificity::ProteinNTerm | TermSpecificity::ProteinCTerm
                    )
                {
                    m.full_id().into()
                } else {
                    m.name().into()
                };
                format!("({name})")
            }
        };
        let mut result = String::new();
        if let Some(m) = &self.n_terminal {
            result.push('.');
            result.push_str(&annotation(m));
        }
        for (residue, modification) in self.sequence.chars().zip(&self.residue_modifications) {
            result.push(residue);
            if let Some(m) = modification {
                result.push_str(&annotation(m));
            }
        }
        if let Some(m) = &self.c_terminal {
            result.push('.');
            result.push_str(&annotation(m));
        }
        result
    }
    /// b/y fragments at internal bonds, charges 1..=max_charge. Mass-only
    /// annotations are supported; unresolved residue masses return Unsupported.
    pub fn fragment_ions(&self, max_charge: u8) -> Result<Vec<FragmentIon>> {
        if max_charge == 0 {
            return Err(invalid("fragment charge must be positive"));
        }
        if self.len() < 2 {
            return Ok(Vec::new());
        }
        self.mono_mass()?;
        let capacity = (self.len() - 1)
            .checked_mul(usize::from(max_charge))
            .and_then(|n| n.checked_mul(2))
            .ok_or_else(|| invalid("fragment count overflow"))?;
        let mut ions = Vec::new();
        ions.try_reserve(capacity)
            .map_err(|_| invalid("fragment allocation exceeds capacity"))?;
        // Accumulate suffixes independently. Subtracting a huge prefix from
        // the full peptide would erase a small, representable complementary ion.
        let mut suffix_masses = Vec::new();
        suffix_masses
            .try_reserve(self.len())
            .map_err(|_| invalid("fragment allocation exceeds capacity"))?;
        let mut suffix_mass = water().mono_mass()
            + self
                .c_terminal
                .as_ref()
                .map(SequenceModification::diff_mono_mass)
                .transpose()?
                .unwrap_or(0.0);
        for index in (0..self.len()).rev() {
            suffix_mass += self
                .internal_mass(index)?
                .ok_or_else(|| unknown("fragment contains an unresolved residue mass"))?;
            suffix_masses.push(suffix_mass);
        }
        let mut b_mass = self
            .n_terminal
            .as_ref()
            .map(SequenceModification::diff_mono_mass)
            .transpose()?
            .unwrap_or(0.0);
        for offset in 0..self.len() - 1 {
            b_mass += self
                .internal_mass(offset)?
                .ok_or_else(|| unknown("fragment contains an unresolved residue mass"))?;
            valid_mass(b_mass)?;
            let y_mass = valid_mass(suffix_masses[self.len() - offset - 2])?;
            let cut = offset + 1;
            for charge in 1..=max_charge {
                let protonation = f64::from(charge) * PROTON_MASS_U;
                ions.push(FragmentIon {
                    series: IonSeries::B,
                    ordinal: cut,
                    charge,
                    mz: (b_mass + protonation) / f64::from(charge),
                });
                ions.push(FragmentIon {
                    series: IonSeries::Y,
                    ordinal: self.len() - cut,
                    charge,
                    mz: (y_mass + protonation) / f64::from(charge),
                });
            }
        }
        Ok(ions)
    }
    fn internal_formula(&self, index: usize) -> Result<Option<EmpiricalFormula>> {
        let base = residue_composition(self.sequence.as_bytes()[index]).map(composition_formula);
        match &self.residue_modifications[index] {
            Some(SequenceModification::MassTag(_)) => Ok(None),
            Some(SequenceModification::Known(m)) if !m.diff_formula().is_empty() => base
                .map(|base| base.checked_add(m.diff_formula()))
                .transpose(),
            Some(SequenceModification::Known(m)) if missing_modification_formula(m) => {
                // Absolute formulas describe a free residue. They can restore a
                // known composition even when the original B/Z/X was unresolved.
                m.absolute_formula()
                    .map(|formula| formula.checked_sub(&water()))
                    .transpose()
            }
            _ => Ok(base),
        }
    }
    pub(crate) fn internal_mass(&self, index: usize) -> Result<Option<f64>> {
        let composition = residue_composition(self.sequence.as_bytes()[index]);
        let base = composition.map(|c| composition_formula(c).mono_mass());
        let mass = match &self.residue_modifications[index] {
            Some(SequenceModification::MassTag(tag)) => tag.residue_mono_mass,
            Some(SequenceModification::Known(m)) if missing_modification_formula(m) => {
                // Source absolute values are free-residue masses. A nonzero
                // value overrides even an unknown base residue; otherwise the
                // declared delta is added before removing full-residue water.
                let water = water();
                let full = if let Some(formula) = m.absolute_formula() {
                    Some(formula.mono_mass())
                } else if m.mono_mass() != 0.0 {
                    Some(m.mono_mass())
                } else {
                    composition
                        .map(|c| {
                            Ok::<_, Error>(
                                composition_formula(c).checked_add(&water)?.mono_mass()
                                    + m.diff_mono_mass(),
                            )
                        })
                        .transpose()?
                };
                full.map(|mass| mass - water.mono_mass())
            }
            Some(SequenceModification::Known(m)) => {
                // A delta formula takes precedence over stored absolute masses.
                // Empty formula plus zero differences is the source no-op case.
                base.map(|base| {
                    if m.diff_formula().is_empty() {
                        base
                    } else {
                        base + m.diff_formula().mono_mass()
                    }
                })
            }
            None => base,
        };
        mass.map(valid_mass).transpose()
    }
    fn rebuild_chemistry(&mut self) -> Result<()> {
        let mut formula = if self.is_empty() {
            Some(EmpiricalFormula::default())
        } else {
            Some(water())
        };
        // Combine terminal deltas before adding residue masses. This retains
        // the peptide mass when large, opposite terminal annotations cancel.
        let mut terminal_delta = 0.0;
        for modification in [&self.n_terminal, &self.c_terminal].into_iter().flatten() {
            terminal_delta += modification.diff_mono_mass()?;
        }
        if !terminal_delta.is_finite() {
            return Err(invalid("terminal annotation mass sum overflows"));
        }
        let mut mass = Some(
            if self.is_empty() {
                0.0
            } else {
                water().mono_mass()
            } + terminal_delta,
        );
        for index in 0..self.len() {
            formula = match (formula, self.internal_formula(index)?) {
                (Some(formula), Some(internal)) => Some(formula.checked_add(&internal)?),
                _ => None,
            };
            let internal = self.internal_mass(index)?;
            mass = match (mass, internal) {
                (Some(a), Some(b)) => {
                    let sum = a + b;
                    if !sum.is_finite() {
                        return Err(invalid("sequence mass sum overflows"));
                    }
                    Some(sum)
                }
                _ => None,
            };
        }
        for modification in [&self.n_terminal, &self.c_terminal].into_iter().flatten() {
            formula = match (formula, modification.known()) {
                (Some(formula), Some(m)) if !missing_modification_formula(m) => {
                    Some(formula.checked_add(m.diff_formula())?)
                }
                _ => None,
            };
        }
        if let Some(formula) = &formula {
            if formula.atoms.values().any(|&count| count < 0) {
                return Err(invalid(
                    "modifications produce negative peptide atom counts",
                ));
            }
        }
        self.formula = formula;
        self.mono_mass = mass.map(valid_mass).transpose()?;
        Ok(())
    }
}

// Residue::setModification only changes a residue if a difference is specified.
// Absolute-only zero-difference provider records describe the existing residue.
fn missing_modification_formula(modification: &ResidueModification) -> bool {
    modification.diff_formula().is_empty()
        && (modification.diff_mono_mass() != 0.0 || modification.diff_average_mass() != 0.0)
}

// Conservative BTreeMap comparison/node and scratch accounting, matching the
// existing RNA formula boundary. A sparse map still needs a complete root node.
fn reserve_formula(entries: usize, work: &mut usize, bytes: &mut usize) -> Result<()> {
    let depth = (usize::BITS - entries.max(1).leading_zeros()) as usize + 1;
    *work = work
        .checked_sub(entries.saturating_mul(depth).saturating_mul(12))
        .ok_or_else(|| invalid("peptide formula work limit exceeded"))?;
    *bytes = bytes
        .checked_sub(entries.saturating_mul(128).saturating_add(512))
        .ok_or_else(|| invalid("peptide formula allocation limit exceeded"))?;
    Ok(())
}

fn consume_mass_work(amount: usize, work: &mut usize) -> Result<()> {
    *work = work
        .checked_sub(amount)
        .ok_or_else(|| invalid("peptide mass work limit exceeded"))?;
    Ok(())
}

fn finite_fragment_mass(mass: f64) -> Result<f64> {
    if mass.is_finite() {
        Ok(mass)
    } else {
        Err(invalid("fragment mass sum overflows"))
    }
}

fn add_formula(
    result: &mut EmpiricalFormula,
    other: &EmpiricalFormula,
    work: &mut usize,
    bytes: &mut usize,
) -> Result<()> {
    reserve_formula(
        result.atoms.len().saturating_add(other.atoms.len()),
        work,
        bytes,
    )?;
    *result = result.checked_add(other)?;
    Ok(())
}

pub(crate) const GENERATION_FORMULA_ENTRY_BYTES: usize =
    std::mem::size_of::<(super::Atom, i32)>() + 3 * std::mem::size_of::<usize>();

#[derive(Clone, Copy)]
enum Attachment {
    NTerm,
    CTerm,
    Residue(usize),
}
#[derive(Clone, Copy)]
enum Annotation<'a> {
    Named(&'a str),
    Mass(&'a str),
}

impl FromStr for AASequence {
    type Err = Error;
    fn from_str(input: &str) -> Result<Self> {
        Self::parse_with_registry(input, ModificationsDB::global())
    }
}
impl AASequence {
    /// Coverage parsing shares the graph's allowance. Preflight conservatively
    /// bounds registry lookups, temporary candidates and cumulative formula copies.
    pub(crate) fn parse_with_budget(
        input: &str,
        db: &ModificationsDB,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<Self> {
        let consume = |remaining: &mut usize, amount: usize| -> Result<()> {
            *remaining = remaining
                .checked_sub(amount)
                .ok_or_else(|| invalid("peptide parent parsing resource limit exceeded"))?;
            Ok(())
        };
        consume(remaining_work, input.len().saturating_add(1))?;
        let annotations = input.bytes().filter(|b| matches!(b, b'(' | b'[')).count();
        let (mut max_id, mut max_atoms, mut keys) = (0usize, 0usize, db.len());
        if annotations != 0 {
            consume(remaining_work, db.len())?;
            for record in db.entries() {
                max_id = max_id.max(record.full_id().len());
                max_atoms = max_atoms.max(
                    record.diff_formula().stored_atom_types().saturating_add(
                        record
                            .absolute_formula()
                            .map_or(0, |f| f.stored_atom_types()),
                    ),
                );
                keys = keys.saturating_add(record.synonyms().len().saturating_add(8));
            }
        }
        let nodes = 6usize.saturating_add(annotations.saturating_mul(max_atoms));
        let copies = input.len().saturating_add(1).saturating_mul(nodes);
        let lookup = db
            .len()
            .saturating_mul(max_id.saturating_add(16))
            .saturating_mul(16)
            .saturating_add(
                input
                    .len()
                    .saturating_add(1)
                    .saturating_mul(keys.max(1).ilog2() as usize + 1)
                    .saturating_mul(288),
            );
        consume(
            remaining_work,
            copies
                .saturating_mul(64)
                .saturating_add(annotations.saturating_mul(lookup)),
        )?;
        consume(
            remaining_bytes,
            input
                .len()
                .saturating_mul(512)
                .saturating_add(copies.saturating_mul(128))
                .saturating_add(annotations.saturating_mul(db.len()).saturating_mul(128)),
        )?;
        Self::parse_with_registry(input, db)
    }

    /// Parse using a supplied registry. Resolved records are shared owned handles,
    /// so neither the registry nor input text must outlive the returned peptide.
    pub fn parse_with_registry(input: &str, db: &ModificationsDB) -> Result<Self> {
        let mut result = Self::default();
        let mut pending = Vec::new();
        let bytes = input.as_bytes();
        let mut offset = 0;
        let mut terminal = None;
        let mut c_terminal_started = false;
        while offset < bytes.len() {
            match bytes[offset] {
                b'.' => {
                    terminal = Some(if result.sequence.is_empty() {
                        Attachment::NTerm
                    } else {
                        c_terminal_started = true;
                        Attachment::CTerm
                    });
                    offset += 1;
                }
                b'n' if result.sequence.is_empty() && offset == 0 => {
                    terminal = Some(Attachment::NTerm);
                    offset += 1;
                }
                b'c' if !result.sequence.is_empty() => {
                    terminal = Some(Attachment::CTerm);
                    c_terminal_started = true;
                    offset += 1;
                }
                bracket @ (b'(' | b'[') => {
                    let start = offset + 1;
                    let mut depth = 1;
                    offset += 1;
                    while offset < bytes.len() && depth > 0 {
                        if bracket == b'(' {
                            match bytes[offset] {
                                b'(' => depth += 1,
                                b')' => depth -= 1,
                                _ => {}
                            }
                        } else if bytes[offset] == b']' {
                            depth = 0;
                        }
                        if depth > 128 {
                            return Err(invalid("modification parentheses nested too deeply"));
                        }
                        if depth > 0 {
                            offset += 1;
                        }
                    }
                    if depth != 0 {
                        return Err(invalid("unterminated sequence annotation"));
                    }
                    let text = &input[start..offset];
                    if text.is_empty() {
                        return Err(invalid("empty sequence annotation"));
                    }
                    let attachment = terminal.take().unwrap_or_else(|| {
                        if result.sequence.is_empty() {
                            Attachment::NTerm
                        } else {
                            Attachment::Residue(result.len() - 1)
                        }
                    });
                    let annotation = if bracket == b'(' {
                        Annotation::Named(text)
                    } else {
                        Annotation::Mass(text)
                    };
                    pending.push((attachment, annotation));
                    offset += 1;
                }
                residue
                    if residue_composition(residue).is_some()
                        || matches!(residue, b'B' | b'Z' | b'X') =>
                {
                    if c_terminal_started {
                        return Err(invalid("residue after C-terminal delimiter"));
                    }
                    terminal = None;
                    result.sequence.push(char::from(residue));
                    result.residue_modifications.push(None);
                    offset += 1;
                }
                _ => {
                    return Err(invalid(format!(
                        "unsupported amino-acid syntax at byte {offset}"
                    )));
                }
            }
        }
        if result.is_empty() && (!pending.is_empty() || !input.is_empty()) {
            return Err(invalid("terminal notation requires a peptide"));
        }
        for (attachment, annotation) in pending {
            let (destination, modification) = match attachment {
                Attachment::NTerm | Attachment::CTerm => {
                    let n_terminal = matches!(attachment, Attachment::NTerm);
                    let residue = result.residue(if n_terminal { 0 } else { result.len() - 1 })?;
                    let modification = match annotation {
                        Annotation::Named(name) => SequenceModification::Known(resolve_terminal(
                            name, residue, n_terminal, db,
                        )?),
                        Annotation::Mass(text) => resolve_mass(
                            text,
                            residue,
                            if n_terminal {
                                TermSpecificity::NTerm
                            } else {
                                TermSpecificity::CTerm
                            },
                            db,
                        )?,
                    };
                    (attachment, modification)
                }
                Attachment::Residue(index) => {
                    let residue = result.residue(index)?;
                    match annotation {
                        Annotation::Named(name) => {
                            let anywhere =
                                db.find(name, Some(residue), Some(TermSpecificity::Anywhere));
                            if !anywhere.is_empty() {
                                (
                                    attachment,
                                    SequenceModification::Known(db.get_modification_handle(
                                        name,
                                        Some(residue),
                                        Some(TermSpecificity::Anywhere),
                                    )?),
                                )
                            } else if index == 0
                                && result.n_terminal.is_none()
                                && resolve_terminal(name, residue, true, db).is_ok()
                            {
                                (
                                    Attachment::NTerm,
                                    SequenceModification::Known(resolve_terminal(
                                        name, residue, true, db,
                                    )?),
                                )
                            } else if index + 1 == result.len() && result.c_terminal.is_none() {
                                (
                                    Attachment::CTerm,
                                    SequenceModification::Known(resolve_terminal(
                                        name, residue, false, db,
                                    )?),
                                )
                            } else {
                                return Err(invalid(format!(
                                    "modification {name:?} does not apply to residue {index}"
                                )));
                            }
                        }
                        Annotation::Mass(text) => (
                            attachment,
                            resolve_mass(text, residue, TermSpecificity::Anywhere, db)?,
                        ),
                    }
                }
            };
            let slot = result.slot(destination);
            if slot.is_some() {
                return Err(invalid("duplicate annotation at a residue or terminus"));
            }
            *slot = Some(modification);
        }
        result.rebuild_chemistry()?;
        Ok(result)
    }
}

fn parse_mass(text: &str) -> Result<(f64, bool)> {
    // Scientific notation would change the source precision calculation (which
    // counts characters after '.') and is deliberately not accepted as a tag.
    if text.is_empty() || text.len() > 1024 {
        return Err(invalid("mass tag must contain 1..=1024 characters"));
    }
    let delta = matches!(text.as_bytes()[0], b'+' | b'-');
    let number = if delta { &text[1..] } else { text };
    if number.is_empty()
        || number.bytes().filter(|b| *b == b'.').count() > 1
        || !number.bytes().all(|b| b.is_ascii_digit() || b == b'.')
        || !number.bytes().any(|b| b.is_ascii_digit())
    {
        return Err(invalid(
            "mass tag requires a finite decimal number without whitespace or exponent",
        ));
    }
    let mass = text
        .parse::<f64>()
        .map_err(|_| invalid("invalid mass tag number"))?;
    if !mass.is_finite() || (mass == 0.0 && number.bytes().any(|b| b.is_ascii_digit() && b != b'0'))
    {
        return Err(invalid(
            "mass tag number must be finite and representable without underflow",
        ));
    }
    Ok((mass, delta))
}
fn lookup_mass(
    text: &str,
    delta: f64,
    residue: char,
    term: TermSpecificity,
    db: &ModificationsDB,
) -> Result<Option<Arc<ResidueModification>>> {
    if let Some(point) = text.find('.') {
        let decimals = i32::try_from(text.len() - point - 1)
            .map_err(|_| invalid("mass tag precision overflows"))?;
        db.best_by_mass_handle(delta, 10_f64.powi(-decimals), Some(residue), Some(term))
    } else {
        Ok(db
            .search_by_mass_handles(delta, 0.5, Some(residue), Some(term))?
            .first()
            .cloned())
    }
}
fn resolve_mass(
    text: &str,
    residue: char,
    term: TermSpecificity,
    db: &ModificationsDB,
) -> Result<SequenceModification> {
    let (mass, is_delta) = parse_mass(text)?;
    let anywhere = term == TermSpecificity::Anywhere;
    let base = if anywhere {
        residue_composition(residue as u8).map(|c| composition_formula(c).mono_mass())
    } else {
        Some(terminal_base(term))
    };
    if is_delta && base.is_none() {
        return Err(invalid(
            "a mass difference cannot modify an unresolved residue mass",
        ));
    }
    let delta = if is_delta {
        Some(mass)
    } else {
        base.map(|base| mass - base)
    };
    // B/Z/X placeholders are not chemical masses. Preserve explicit absolute
    // tags rather than looking up a delta against a fictitious zero mass.
    if let Some(delta) = delta {
        if !delta.is_finite() {
            return Err(invalid("mass tag delta overflow"));
        }
        if let Some(known) = lookup_mass(text, delta, residue, term, db)? {
            return Ok(SequenceModification::Known(known));
        }
    }
    let internal = if anywhere {
        Some(valid_mass(if is_delta {
            base.unwrap() + mass
        } else {
            mass
        })?)
    } else {
        None
    };
    let full_id = if anywhere {
        format!("{residue}[{text}]")
    } else if matches!(term, TermSpecificity::NTerm | TermSpecificity::ProteinNTerm) {
        format!(".n[{text}]")
    } else {
        format!(".c[{text}]")
    };
    Ok(SequenceModification::MassTag(MassTag {
        input: text.into(),
        full_id: full_id.into(),
        mass,
        is_delta,
        delta_mono_mass: delta,
        residue_mono_mass: internal,
        origin: if anywhere { Some(residue) } else { None },
        term,
    }))
}
fn resolve_terminal(
    name: &str,
    residue: char,
    n_terminal: bool,
    db: &ModificationsDB,
) -> Result<Arc<ResidueModification>> {
    let types = if n_terminal {
        [TermSpecificity::NTerm, TermSpecificity::ProteinNTerm]
    } else {
        [TermSpecificity::CTerm, TermSpecificity::ProteinCTerm]
    };
    for term in types {
        if !db.find(name, Some(residue), Some(term)).is_empty() {
            return db.get_modification_handle(name, Some(residue), Some(term));
        }
    }
    Err(invalid(format!(
        "no modification {name:?} matches the requested terminus"
    )))
}
impl fmt::Display for AASequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.annotated_string(false))
    }
}

#[cfg(test)]
mod parent_parse_budget_tests {
    use super::*;

    #[test]
    fn repeated_parent_parses_share_work_and_bytes_without_changing_chemistry() {
        let db = ModificationsDB::global();
        const WORK: usize = 50_000_000;
        const BYTES: usize = 256 * 1024 * 1024;
        for text in ["PEPTIDE", "M(Oxidation)PEPTIDE", "A[+12.345678]G"] {
            let expected = AASequence::parse_with_registry(text, db).unwrap();
            let (mut work, mut bytes) = (WORK, BYTES);
            assert_eq!(
                AASequence::parse_with_budget(text, db, &mut work, &mut bytes).unwrap(),
                expected
            );
            let (used_work, used_bytes) = (WORK - work, BYTES - bytes);
            assert!(used_work > 0 && used_bytes > 0);
            for constrain_work in [true, false] {
                let (mut work, mut bytes) = if constrain_work {
                    (used_work * 2 - 1, BYTES)
                } else {
                    (WORK, used_bytes * 2 - 1)
                };
                assert_eq!(
                    AASequence::parse_with_budget(text, db, &mut work, &mut bytes).unwrap(),
                    expected
                );
                assert!(AASequence::parse_with_budget(text, db, &mut work, &mut bytes).is_err());
            }
        }
    }
}

#[cfg(test)]
mod fragment_formula_budget_tests {
    use super::*;

    #[test]
    fn formula_calls_share_work_and_allocation_counters() {
        let peptide = AASequence::parse("(Acetyl)AC(Carbamidomethyl)M(Oxidation)K").unwrap();
        let expected = peptide.formula_for(PeptideFragmentType::Full, 2).unwrap();
        let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
        assert_eq!(
            peptide
                .formula_for_with_budget(PeptideFragmentType::Full, 2, &mut work, &mut bytes)
                .unwrap(),
            expected
        );
        let (used_work, used_bytes) = (50_000_000 - work, 256 * 1024 * 1024 - bytes);
        assert!(used_work > 0 && used_bytes > 0);
        for constrain_work in [true, false] {
            let (mut work, mut bytes) = if constrain_work {
                (2 * used_work - 1, 256 * 1024 * 1024)
            } else {
                (50_000_000, 2 * used_bytes - 1)
            };
            assert_eq!(
                peptide
                    .formula_for_with_budget(PeptideFragmentType::Full, 2, &mut work, &mut bytes)
                    .unwrap(),
                expected
            );
            let before = (work, bytes);
            assert!(
                peptide
                    .formula_for_with_budget(PeptideFragmentType::Full, 2, &mut work, &mut bytes)
                    .is_err()
            );
            assert!(work < before.0 || bytes < before.1);
            assert_eq!(
                peptide.formula_for(PeptideFragmentType::Full, 2).unwrap(),
                expected
            );
        }
    }

    #[test]
    fn unknown_composition_and_early_budget_failure_leave_sequence_unchanged() {
        let peptide = AASequence::parse("AX").unwrap();
        let original = peptide.clone();
        let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
        assert!(matches!(
            peptide.formula_for_with_budget(PeptideFragmentType::Full, 0, &mut work, &mut bytes),
            Err(Error::Unsupported(_))
        ));
        assert!(work < 50_000_000 && bytes < 256 * 1024 * 1024);
        assert!(
            peptide
                .formula_for_with_budget(PeptideFragmentType::Full, 0, &mut 0, &mut 0)
                .is_err()
        );
        assert_eq!(peptide, original);
    }

    #[test]
    fn empty_typed_terminal_state_needs_no_formula_budget() {
        let modification = ModificationsDB::global()
            .get_modification_handle("Acetyl", None, Some(TermSpecificity::NTerm))
            .unwrap();
        let empty = AASequence::default()
            .with_resolved_modifications(&[], Some(modification), None)
            .unwrap();
        assert!(empty.n_terminal_modification().is_some());
        let (mut work, mut bytes) = (0, 0);
        assert_eq!(
            empty
                .formula_for_with_budget(PeptideFragmentType::Full, i32::MAX, &mut work, &mut bytes)
                .unwrap(),
            EmpiricalFormula::default()
        );
        assert_eq!((work, bytes), (0, 0));
        assert_eq!(
            empty
                .mono_mass_for_with_budget(
                    PeptideFragmentType::Full,
                    i32::MAX,
                    &mut work,
                    &mut bytes
                )
                .unwrap(),
            0.0
        );
        assert_eq!((work, bytes), (0, 0));
    }

    #[test]
    fn scalar_fragment_masses_share_limits_including_mass_only_records() {
        let db = ModificationsDB::from_records(vec![
            ResidueModification::from_record(super::super::ModificationRecord {
                name: "Delta".into(),
                origin: Some('M'),
                diff_mono_mass: 12.5,
                ..Default::default()
            })
            .unwrap(),
        ])
        .unwrap();
        for peptide in [
            AASequence::parse("AG").unwrap(),
            AASequence::parse("AM(Oxidation)").unwrap(),
            AASequence::parse("AX[999]").unwrap(),
            AASequence::parse_with_registry("AM(Delta)", &db).unwrap(),
        ] {
            let (mut work, mut bytes) = (50_000_000, 256 * 1024 * 1024);
            let mass = peptide
                .mono_mass_for_with_budget(PeptideFragmentType::Full, 2, &mut work, &mut bytes)
                .unwrap();
            let (used_work, used_bytes) = (50_000_000 - work, 256 * 1024 * 1024 - bytes);
            for constrain_work in [true, false] {
                let (mut work, mut bytes) = if constrain_work {
                    (used_work * 2 - 1, 256 * 1024 * 1024)
                } else {
                    (50_000_000, used_bytes * 2 - 1)
                };
                assert_eq!(
                    peptide
                        .mono_mass_for_with_budget(
                            PeptideFragmentType::Full,
                            2,
                            &mut work,
                            &mut bytes
                        )
                        .unwrap(),
                    mass
                );
                assert!(
                    peptide
                        .mono_mass_for_with_budget(
                            PeptideFragmentType::Full,
                            2,
                            &mut work,
                            &mut bytes
                        )
                        .is_err()
                );
            }
        }
    }
}
