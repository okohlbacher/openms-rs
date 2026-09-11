// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::resolution::{Budget, Resolver};
use super::*;
use crate::chemistry::{
    EmpiricalFormula, ModificationsDB, PROTON_MASS_U, composition_formula, residue_composition,
};
use std::{collections::BTreeSet, mem::size_of};

pub const MAX_PROFORMA_MASS_WORK: usize = MAX_PROFORMA_RESOLUTION_WORK;
pub const MAX_PROFORMA_MASS_BYTES: usize = MAX_PROFORMA_RESOLUTION_BYTES;
pub const MAX_PROFORMA_MASS_ITEMS: usize = MAX_PROFORMA_RESOLUTION_ITEMS;
pub const MAX_PROFORMA_MASS_TEXT_BYTES: usize = MAX_PROFORMA_RESOLUTION_TEXT_BYTES;

/// A source mass/issue/predicate value together with its owned resolution warnings.
#[derive(Clone, Debug, PartialEq)]
pub struct MassEvaluation<T> {
    pub value: T,
    pub warnings: Vec<ResolutionWarning>,
}

/// Completed optional mass/mz calculation. `None` denotes scientific unavailability;
/// invalid numerical state or a resource failure is a separate `Err`.
#[derive(Clone, Debug, PartialEq)]
pub struct MassAttempt {
    pub value: Option<f64>,
    pub issues: Vec<ConversionIssue>,
    pub warnings: Vec<ResolutionWarning>,
}

impl Peptidoform {
    pub fn mass_calculation_issues(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<MassEvaluation<Vec<ConversionIssue>>> {
        issues(Input::Chain(self), registry)
    }
    pub fn can_calculate_mass(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<MassEvaluation<bool>> {
        let result = self.mass_calculation_issues(registry)?;
        Ok(MassEvaluation {
            value: result.value.is_empty(),
            warnings: result.warnings,
        })
    }
    /// Source ordered neutral mass, including its documented range/ambiguity and
    /// crosslink omissions. Errors leave the registry and AST unchanged.
    pub fn mono_mass(&self, registry: &mut ModificationsDB) -> Result<MassEvaluation<f64>> {
        scalar(Input::Chain(self), registry, None)
    }
    pub fn mz(&self, charge: i32, registry: &mut ModificationsDB) -> Result<MassEvaluation<f64>> {
        scalar(Input::Chain(self), registry, Some(Charge::Explicit(charge)))
    }
    pub fn try_mono_mass(&self, registry: &mut ModificationsDB) -> Result<MassAttempt> {
        run(
            Input::Chain(self),
            registry,
            Mode::Try,
            None,
            &mut Budget::default(),
        )
    }
    pub fn try_mz(&self, charge: i32, registry: &mut ModificationsDB) -> Result<MassAttempt> {
        run(
            Input::Chain(self),
            registry,
            Mode::Try,
            Some(Charge::Explicit(charge)),
            &mut Budget::default(),
        )
    }
}
impl PeptidoformIon {
    pub fn mass_calculation_issues(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<MassEvaluation<Vec<ConversionIssue>>> {
        issues(Input::Ion(self), registry)
    }
    pub fn can_calculate_mass(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<MassEvaluation<bool>> {
        let result = self.mass_calculation_issues(registry)?;
        Ok(MassEvaluation {
            value: result.value.is_empty(),
            warnings: result.warnings,
        })
    }
    pub fn mono_mass(&self, registry: &mut ModificationsDB) -> Result<MassEvaluation<f64>> {
        scalar(Input::Ion(self), registry, None)
    }
    pub fn mz(&self, registry: &mut ModificationsDB) -> Result<MassEvaluation<f64>> {
        scalar(
            Input::Ion(self),
            registry,
            Some(Charge::Stored(self.charge.as_ref())),
        )
    }
    pub fn try_mono_mass(&self, registry: &mut ModificationsDB) -> Result<MassAttempt> {
        run(
            Input::Ion(self),
            registry,
            Mode::Try,
            None,
            &mut Budget::default(),
        )
    }
    pub fn try_mz(&self, registry: &mut ModificationsDB) -> Result<MassAttempt> {
        run(
            Input::Ion(self),
            registry,
            Mode::Try,
            Some(Charge::Stored(self.charge.as_ref())),
            &mut Budget::default(),
        )
    }
}

#[derive(Clone, Copy)]
enum Input<'a> {
    Chain(&'a Peptidoform),
    Ion(&'a PeptidoformIon),
}
impl<'a> Input<'a> {
    fn chains(self) -> &'a [Peptidoform] {
        match self {
            Self::Chain(chain) => std::slice::from_ref(chain),
            Self::Ion(ion) => &ion.chains,
        }
    }
    fn chimeric(self) -> bool {
        matches!(self, Self::Ion(ion) if ion.is_chimeric)
    }
}
#[derive(Clone, Copy, PartialEq)]
enum Mode {
    Issues,
    Strict,
    Try,
}
#[derive(Clone, Copy)]
enum Charge<'a> {
    Explicit(i32),
    Stored(Option<&'a ChargeState>),
}

fn issues(
    input: Input<'_>,
    registry: &mut ModificationsDB,
) -> Result<MassEvaluation<Vec<ConversionIssue>>> {
    let report = run(input, registry, Mode::Issues, None, &mut Budget::default())?;
    Ok(MassEvaluation {
        value: report.issues,
        warnings: report.warnings,
    })
}
fn scalar(
    input: Input<'_>,
    registry: &mut ModificationsDB,
    charge: Option<Charge<'_>>,
) -> Result<MassEvaluation<f64>> {
    let report = run(
        input,
        registry,
        Mode::Strict,
        charge,
        &mut Budget::default(),
    )?;
    Ok(MassEvaluation {
        value: report.value.expect("strict mass result"),
        warnings: report.warnings,
    })
}

fn run(
    input: Input<'_>,
    registry: &mut ModificationsDB,
    mode: Mode,
    charge: Option<Charge<'_>>,
    budget: &mut Budget,
) -> Result<MassAttempt> {
    let mut report = MassAttempt {
        value: None,
        issues: Vec::new(),
        warnings: Vec::new(),
    };
    // Source charge diagnostics precede all sequence work, including empty ions.
    let charge = match charge {
        Some(charge) => match charge_value(charge, budget)? {
            Some(value) if value != 0 => Some(value),
            value => {
                let message = if value.is_none() {
                    "No charge state specified"
                } else {
                    "Charge state is zero"
                };
                if mode == Mode::Strict {
                    return Err(invalid(if value.is_none() {
                        "Cannot calculate m/z: no charge state specified"
                    } else {
                        "Cannot calculate m/z: charge state is zero"
                    }));
                }
                issue(
                    &mut report.issues,
                    ConversionIssueType::UnsupportedFeature,
                    message,
                    Some(0),
                    budget,
                )?;
                return Ok(report);
            }
        },
        None => None,
    };
    let chains = input.chains();
    budget.items(chains.len())?;
    if mode == Mode::Try && matches!(input, Input::Ion(_)) {
        if chains.is_empty() {
            report.value = Some(finish_mz(0., charge)?);
            return Ok(report);
        }
        if input.chimeric() {
            issue(
                &mut report.issues,
                ConversionIssueType::UnsupportedFeature,
                "Cannot calculate single mass for chimeric spectra.",
                Some(0),
                budget,
            )?;
            return Ok(report);
        }
    }
    let mut resolver = Resolver::new(registry, budget);
    let mut masses = ResidueMasses::new(resolver.budget)?;
    let mut labels = BTreeSet::new();
    if mode == Mode::Try && matches!(input, Input::Chain(_)) {
        // CPP-020: keep source copy A for calculation and resolve copy B
        // inside check_chain. Formula interning can make the passes differ.
        let mut first = clone_chain(&chains[0], resolver.budget)?;
        resolver.resolve_chain(&mut first)?;
        check_chain(&first, None, &mut report.issues, &mut resolver, &mut masses)?;
        if report.issues.is_empty() {
            report.value = Some(finish_mz(
                calculate(&first, &mut labels, resolver.budget, &mut masses)?,
                charge,
            )?);
        }
    } else {
        for (index, chain) in chains.iter().enumerate() {
            let prefix = matches!(input, Input::Ion(_)).then_some(index);
            if mode == Mode::Try {
                let mut first = clone_chain(chain, resolver.budget)?;
                resolver.resolve_chain(&mut first)?;
                check_chain(
                    &first,
                    prefix,
                    &mut report.issues,
                    &mut resolver,
                    &mut masses,
                )?;
            } else {
                check_chain(
                    chain,
                    prefix,
                    &mut report.issues,
                    &mut resolver,
                    &mut masses,
                )?;
            }
        }
        if mode == Mode::Strict && !report.issues.is_empty() {
            resolver
                .budget
                .allocate(report.issues[0].description.len().saturating_add(32))?;
            resolver
                .budget
                .consume(report.issues[0].description.len().saturating_add(32))?;
            return Err(Error::InvalidValue(format!(
                "Cannot calculate mass: {}",
                report.issues[0].description
            )));
        }
        if mode != Mode::Issues && report.issues.is_empty() {
            // Strict source checks every chain's issues before chimeric state,
            // and returns zero for no chains even if is_chimeric is set.
            if !chains.is_empty() && input.chimeric() {
                return Err(invalid(
                    "Cannot calculate single mass for chimeric spectra.",
                ));
            }
            let mut total = 0.;
            for chain in chains {
                let mut resolved = clone_chain(chain, resolver.budget)?;
                resolver.resolve_chain(&mut resolved)?;
                total = finite(
                    total + calculate(&resolved, &mut labels, resolver.budget, &mut masses)?,
                )?;
            }
            report.value = Some(finish_mz(total, charge)?);
        }
    }
    let (staged, warnings) = resolver.finish();
    report.warnings = warnings;
    // Scientific-unavailable Ok reports are completed source operations. Any
    // Err above instead drops staged additions and leaves the caller untouched.
    if let Some(staged) = staged {
        *registry = staged;
    }
    Ok(report)
}

fn charge_value(charge: Charge<'_>, budget: &mut Budget) -> Result<Option<i32>> {
    match charge {
        Charge::Explicit(value) => Ok(Some(value)),
        Charge::Stored(Some(ChargeState::Simple(value))) => Ok(Some(*value)),
        Charge::Stored(None) => Ok(None),
        Charge::Stored(Some(ChargeState::Adducts(adducts))) => {
            budget.items(adducts.len())?;
            budget.consume(adducts.len().saturating_mul(3))?;
            let mut total = 0_i32;
            for adduct in adducts {
                total = adduct
                    .charge
                    .checked_mul(adduct.occurrence.unwrap_or(1))
                    .and_then(|v| total.checked_add(v))
                    .ok_or_else(|| invalid("ProForma adduct charge overflow"))?;
            }
            Ok(Some(total))
        }
    }
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "ProForma mass calculation requires finite consumed values and results",
        ))
    }
}
fn finish_mz(mass: f64, charge: Option<i32>) -> Result<f64> {
    match charge {
        Some(charge) => finite(
            (finite(mass + f64::from(charge) * PROTON_MASS_U)?) / f64::from(charge.unsigned_abs()),
        ),
        None => finite(mass),
    }
}

struct ResidueMasses {
    cache: [Option<Option<f64>>; 128],
    water: f64,
}
impl ResidueMasses {
    fn new(budget: &mut Budget) -> Result<Self> {
        budget.consume(1024)?;
        budget.allocate(1024)?;
        Ok(Self {
            cache: [None; 128],
            water: composition_formula([0, 2, 0, 1, 0, 0]).mono_mass(),
        })
    }
    fn get(&mut self, amino_acid: char, budget: &mut Budget) -> Result<Option<f64>> {
        budget.consume(1)?;
        if !amino_acid.is_ascii() {
            return Ok(None);
        }
        let index = amino_acid as usize;
        if let Some(value) = self.cache[index] {
            return Ok(value);
        }
        budget.consume(2048)?;
        budget.allocate(2048)?;
        let value = if let Some(composition) = residue_composition(amino_acid as u8) {
            let full = composition_formula(composition)
                .checked_add(&composition_formula([0, 2, 0, 1, 0, 0]))?;
            Some(full.mono_mass() - self.water)
        } else if matches!(amino_acid, 'B' | 'Z' | 'X') {
            Some(-self.water)
        } else {
            None
        };
        self.cache[index] = Some(value);
        Ok(value)
    }
}

fn modification_mass(modification: &Modification, budget: &mut Budget) -> Result<Option<f64>> {
    budget.consume(1)?;
    if let Some(record) = &modification.resolved_mod {
        return Ok(Some(finite(record.diff_mono_mass())?));
    }
    budget.items(modification.alternatives.len())?;
    let Some((first, _)) = modification.alternatives.first() else {
        return Ok(None);
    };
    let tag = modification
        .alternatives
        .iter()
        .find_map(|(tag, _)| {
            (!matches!(
                tag,
                ModificationTag::InfoTag(_) | ModificationTag::PositionConstraint(_)
            ))
            .then_some(tag)
        })
        .unwrap_or(first);
    match tag {
        ModificationTag::MassDelta(delta) => Ok(Some(finite(delta.mass)?)),
        ModificationTag::FormulaTag(formula) => {
            budget.text(&formula.formula_string)?;
            budget.consume(
                formula
                    .formula_string
                    .len()
                    .saturating_add(1)
                    .saturating_mul(256),
            )?;
            budget.allocate(
                formula
                    .formula_string
                    .len()
                    .saturating_add(1)
                    .saturating_mul(512)
                    .saturating_add(2048),
            )?;
            // Unlike resolution, this source helper ignores FormulaTag.charge,
            // accepts an empty formula and retains charge parsed from the text.
            let Ok(formula) = EmpiricalFormula::parse(&formula.formula_string) else {
                return Ok(None);
            };
            budget.consume(formula.stored_atom_types().saturating_mul(32))?;
            Ok(Some(finite(formula.mono_mass())?))
        }
        ModificationTag::InfoTag(_) | ModificationTag::PositionConstraint(_) => Ok(Some(0.)),
        _ => Ok(None),
    }
}
fn issue(
    issues: &mut Vec<ConversionIssue>,
    issue_type: ConversionIssueType,
    description: &str,
    position: Option<usize>,
    budget: &mut Budget,
) -> Result<()> {
    budget.text(description)?;
    budget.allocate(description.len())?;
    budget.reserve(issues)?;
    issues.push(ConversionIssue {
        issue_type,
        description: description.to_owned(),
        position,
    });
    Ok(())
}
fn check_modifications(
    modifications: &[Modification],
    position: Option<usize>,
    issues: &mut Vec<ConversionIssue>,
    budget: &mut Budget,
) -> Result<()> {
    budget.items(modifications.len())?;
    for modification in modifications {
        if modification_mass(modification, budget)?.is_none() {
            budget.consume(128)?;
            budget.allocate(128)?;
            let description = format!(
                "Modification at position {} has no resolvable mass",
                position.unwrap_or(usize::MAX)
            );
            issue(
                issues,
                ConversionIssueType::UnresolvedMod,
                &description,
                position,
                budget,
            )?;
        }
    }
    Ok(())
}
fn unknown(
    amino_acid: char,
    position: usize,
    group: Option<&str>,
    issues: &mut Vec<ConversionIssue>,
    budget: &mut Budget,
) -> Result<()> {
    budget.consume(128)?;
    budget.allocate(128)?;
    let description = match group {
        Some(group) => format!("Unknown amino acid '{amino_acid}' in {group}"),
        None => format!("Unknown amino acid '{amino_acid}' at position {position}"),
    };
    issue(
        issues,
        ConversionIssueType::UnsupportedFeature,
        &description,
        Some(position),
        budget,
    )
}
fn check_chain(
    chain: &Peptidoform,
    prefix: Option<usize>,
    issues: &mut Vec<ConversionIssue>,
    resolver: &mut Resolver<'_>,
    masses: &mut ResidueMasses,
) -> Result<()> {
    let mut resolved = clone_chain(chain, resolver.budget)?;
    resolver.resolve_chain(&mut resolved)?;
    let budget = &mut *resolver.budget;
    let start = issues.len();
    let mut position = 0_usize;
    budget.items(resolved.sequence.len())?;
    for section in &resolved.sequence {
        match section {
            SequenceSection::Element(element) => {
                if masses.get(element.amino_acid, budget)?.is_none() {
                    unknown(element.amino_acid, position, None, issues, budget)?;
                }
                check_modifications(&element.modifications, Some(position), issues, budget)?;
                position += 1;
            }
            SequenceSection::AmbiguousRegion(region) => {
                // CPP-009: source compares base masses and ignores all candidate
                // modifications here, although calculation reads the first's.
                budget.items(region.elements.len())?;
                let mut first = None;
                let mut different = false;
                for element in &region.elements {
                    if let Some(mass) = masses.get(element.amino_acid, budget)? {
                        if let Some(first) = first {
                            if first != mass {
                                different = true;
                            }
                        } else {
                            first = Some(mass);
                        }
                    } else {
                        unknown(
                            element.amino_acid,
                            position,
                            Some("ambiguous region"),
                            issues,
                            budget,
                        )?;
                    }
                }
                if different {
                    issue(
                        issues,
                        ConversionIssueType::AmbiguousRegion,
                        "Ambiguous region contains amino acids with different masses",
                        Some(position),
                        budget,
                    )?;
                }
                position += 1;
            }
            SequenceSection::ModifiedRange(range) => {
                budget.items(range.elements.len())?;
                let beginning = position;
                for element in &range.elements {
                    if masses.get(element.amino_acid, budget)?.is_none() {
                        unknown(element.amino_acid, position, Some("range"), issues, budget)?;
                    }
                    position += 1;
                }
                // CPP-008: inner element modifications are omitted in source.
                check_modifications(&range.modifications, Some(beginning), issues, budget)?;
            }
        }
    }
    check_modifications(&resolved.n_term_mods, Some(0), issues, budget)?;
    check_modifications(
        &resolved.c_term_mods,
        Some(position.saturating_sub(1)),
        issues,
        budget,
    )?;
    budget.items(resolved.unlocalised_mods.len())?;
    for group in &resolved.unlocalised_mods {
        check_modifications(&group.modifications, None, issues, budget)?;
    }
    budget.items(resolved.labile_mods.len())?;
    for labile in &resolved.labile_mods {
        check_modifications(
            std::slice::from_ref(&labile.modification),
            None,
            issues,
            budget,
        )?;
    }
    budget.items(resolved.global_mods.len())?;
    for entry in &resolved.global_mods {
        if let GlobalModEntry::GlobalModification(global) = entry {
            check_modifications(
                std::slice::from_ref(&global.modification),
                None,
                issues,
                budget,
            )?;
        }
    }
    if let Some(index) = prefix {
        budget.consume(issues.len() - start)?;
        for issue in &mut issues[start..] {
            budget.consume(issue.description.len().saturating_add(64))?;
            budget.allocate(issue.description.len().saturating_add(64))?;
            issue.description = format!("Chain {index}: {}", issue.description);
        }
    }
    Ok(())
}

fn add_modification(
    modification: &Modification,
    mass: &mut f64,
    labels: &mut BTreeSet<String>,
    budget: &mut Budget,
) -> Result<()> {
    budget.consume(1)?;
    if let Some((_, Some(label))) = modification.alternatives.first() {
        if label.label_type == LabelType::Crosslink {
            budget.text(&label.identifier)?;
            budget.consume(
                label
                    .identifier
                    .len()
                    .saturating_add(1)
                    .saturating_mul(1024),
            )?;
            if labels.contains(&label.identifier) {
                return Ok(());
            }
            // CPP-015: source reserves the ID before checking whether this
            // endpoint provides chemistry; label-only first endpoints count zero.
            budget.allocate(label.identifier.len().saturating_add(1024))?;
            labels.insert(label.identifier.clone());
        }
    }
    if let Some(value) = modification_mass(modification, budget)? {
        *mass = finite(*mass + value)?;
    }
    Ok(())
}
fn add_group(
    group: &[Modification],
    mass: &mut f64,
    labels: &mut BTreeSet<String>,
    budget: &mut Budget,
) -> Result<()> {
    budget.items(group.len())?;
    for modification in group {
        add_modification(modification, mass, labels, budget)?;
    }
    Ok(())
}
fn calculate(
    chain: &Peptidoform,
    labels: &mut BTreeSet<String>,
    budget: &mut Budget,
    masses: &mut ResidueMasses,
) -> Result<f64> {
    let mut mass = 0.;
    budget.items(chain.sequence.len())?;
    for section in &chain.sequence {
        match section {
            SequenceSection::Element(element) => {
                mass = finite(
                    mass + masses.get(element.amino_acid, budget)?.ok_or_else(|| {
                        invalid("Unknown amino acid during ProForma mass calculation")
                    })?,
                )?;
                add_group(&element.modifications, &mut mass, labels, budget)?;
            }
            SequenceSection::AmbiguousRegion(region) => {
                if let Some(element) = region.elements.first() {
                    mass = finite(
                        mass + masses.get(element.amino_acid, budget)?.ok_or_else(|| {
                            invalid("Unknown ambiguous amino acid during ProForma mass calculation")
                        })?,
                    )?;
                    add_group(&element.modifications, &mut mass, labels, budget)?;
                }
            }
            SequenceSection::ModifiedRange(range) => {
                budget.items(range.elements.len())?;
                for element in &range.elements {
                    mass = finite(
                        mass + masses.get(element.amino_acid, budget)?.ok_or_else(|| {
                            invalid("Unknown range amino acid during ProForma mass calculation")
                        })?,
                    )?;
                }
                add_group(&range.modifications, &mut mass, labels, budget)?;
            }
        }
    }
    mass = finite(mass + masses.water)?;
    add_group(&chain.n_term_mods, &mut mass, labels, budget)?;
    add_group(&chain.c_term_mods, &mut mass, labels, budget)?;
    budget.items(chain.unlocalised_mods.len())?;
    for group in &chain.unlocalised_mods {
        budget.items(group.modifications.len())?;
        for modification in &group.modifications {
            if let Some(value) = modification_mass(modification, budget)? {
                mass = finite(mass + value * f64::from(group.occurrence.unwrap_or(1)))?;
            }
        }
    }
    budget.items(chain.labile_mods.len())?;
    for labile in &chain.labile_mods {
        if let Some(value) = modification_mass(&labile.modification, budget)? {
            mass = finite(mass + value)?;
        }
    }
    budget.items(chain.global_mods.len())?;
    for entry in &chain.global_mods {
        if let GlobalModEntry::GlobalModification(global) = entry {
            if let Some(value) = modification_mass(&global.modification, budget)? {
                let mut count = 0_i32;
                budget.consume(chain.sequence.len())?;
                for section in &chain.sequence {
                    if let SequenceSection::Element(element) = section {
                        for location in &global.locations {
                            budget.consume(1)?;
                            if location.len() == 1
                                && element.amino_acid == char::from(location.as_bytes()[0])
                            {
                                count = count.checked_add(1).ok_or_else(|| {
                                    invalid("ProForma global modification count overflow")
                                })?;
                                break;
                            }
                        }
                    }
                }
                mass = finite(mass + value * f64::from(count))?;
            }
        }
    }
    Ok(mass)
}

// The AST has no recursive variants. Visit its fixed-depth vectors before Clone,
// charging owned strings/slots and clone/drop work; Arc chemistry is shared.
fn copy_text(value: &str, budget: &mut Budget) -> Result<()> {
    budget.text(value)?;
    budget.consume(value.len().saturating_mul(2))?;
    budget.allocate(value.len())
}
fn copy_vector<T>(value: &[T], budget: &mut Budget) -> Result<()> {
    budget.items(value.len())?;
    budget.consume(value.len().saturating_mul(2))?;
    budget.allocate(value.len().saturating_mul(size_of::<T>()))
}
fn copy_modifications(value: &[Modification], budget: &mut Budget) -> Result<()> {
    copy_vector(value, budget)?;
    for modification in value {
        copy_vector(&modification.alternatives, budget)?;
        for (tag, label) in &modification.alternatives {
            match tag {
                ModificationTag::CvAccession(value) => copy_text(&value.accession, budget)?,
                ModificationTag::NamedMod(value) => copy_text(&value.name, budget)?,
                ModificationTag::MassDelta(value) => copy_text(&value.original_text, budget)?,
                ModificationTag::FormulaTag(value) => copy_text(&value.formula_string, budget)?,
                ModificationTag::InfoTag(value) => copy_text(&value.text, budget)?,
                ModificationTag::PositionConstraint(value) => copy_vector(&value.residues, budget)?,
                ModificationTag::GlycanComposition(value) => {
                    copy_vector(&value.components, budget)?;
                    for (component, _) in &value.components {
                        copy_text(
                            match component {
                                GlycanComponent::Name(name) => name,
                                GlycanComponent::Formula(formula) => &formula.formula_string,
                            },
                            budget,
                        )?;
                    }
                }
            }
            if let Some(label) = label {
                copy_text(&label.identifier, budget)?;
            }
        }
    }
    Ok(())
}
fn copy_element(element: &SequenceElement, budget: &mut Budget) -> Result<()> {
    copy_modifications(&element.modifications, budget)
}
fn clone_chain(chain: &Peptidoform, budget: &mut Budget) -> Result<Peptidoform> {
    budget.consume(size_of::<Peptidoform>())?;
    budget.allocate(size_of::<Peptidoform>())?;
    if let Some(name) = &chain.name {
        copy_text(name, budget)?;
    }
    copy_vector(&chain.global_mods, budget)?;
    for entry in &chain.global_mods {
        match entry {
            GlobalModEntry::IsotopeReplacement(value) => copy_text(&value.isotope, budget)?,
            GlobalModEntry::GlobalModification(value) => {
                copy_modifications(std::slice::from_ref(&value.modification), budget)?;
                copy_vector(&value.locations, budget)?;
                for location in &value.locations {
                    copy_text(location, budget)?;
                }
            }
        }
    }
    copy_vector(&chain.unlocalised_mods, budget)?;
    for group in &chain.unlocalised_mods {
        copy_modifications(&group.modifications, budget)?;
    }
    copy_vector(&chain.labile_mods, budget)?;
    for labile in &chain.labile_mods {
        copy_modifications(std::slice::from_ref(&labile.modification), budget)?;
    }
    copy_modifications(&chain.n_term_mods, budget)?;
    copy_modifications(&chain.c_term_mods, budget)?;
    copy_vector(&chain.sequence, budget)?;
    for section in &chain.sequence {
        match section {
            SequenceSection::Element(element) => copy_element(element, budget)?,
            SequenceSection::AmbiguousRegion(region) => {
                copy_vector(&region.elements, budget)?;
                for element in &region.elements {
                    copy_element(element, budget)?;
                }
            }
            SequenceSection::ModifiedRange(range) => {
                copy_vector(&range.elements, budget)?;
                for element in &range.elements {
                    copy_element(element, budget)?;
                }
                copy_modifications(&range.modifications, budget)?;
            }
        }
    }
    if let Some(ChargeState::Adducts(adducts)) = &chain.charge {
        copy_vector(adducts, budget)?;
        for adduct in adducts {
            copy_text(&adduct.formula, budget)?;
        }
    }
    Ok(chain.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn db() -> ModificationsDB {
        ModificationsDB::from_records(Vec::new()).unwrap()
    }
    #[test]
    fn every_final_counter_failure_rolls_back_after_formula_interning() {
        let pf = Peptidoform::parse("M[Formula:Cl101]A[+12]").unwrap();
        let before = pf.clone();
        let mut full = Budget::default();
        let mut original = db();
        assert!(
            run(
                Input::Chain(&pf),
                &mut original,
                Mode::Strict,
                None,
                &mut full
            )
            .is_ok()
        );
        assert_eq!(original.len(), 1);
        let used = [
            MAX_PROFORMA_MASS_WORK - full.work,
            MAX_PROFORMA_MASS_BYTES - full.bytes,
            MAX_PROFORMA_MASS_ITEMS - full.items,
        ];
        for (index, used) in used.into_iter().enumerate() {
            assert!(used > 1);
            let mut limited = Budget::default();
            match index {
                0 => limited.work = used - 1,
                1 => limited.bytes = used - 1,
                _ => limited.items = used - 1,
            }
            let mut original = db();
            assert!(
                run(
                    Input::Chain(&pf),
                    &mut original,
                    Mode::Strict,
                    None,
                    &mut limited
                )
                .is_err()
            );
            assert!(original.is_empty());
            assert_eq!(pf, before);
        }
    }
    #[test]
    fn completed_issue_attempt_commits_but_scalar_error_does_not() {
        let pf = Peptidoform::parse("M[Formula:Cl101]A[missing_name]").unwrap();
        for mode in [Mode::Issues, Mode::Try] {
            let mut original = db();
            let result = run(
                Input::Chain(&pf),
                &mut original,
                mode,
                None,
                &mut Budget::default(),
            )
            .unwrap();
            assert_eq!(original.len(), 1);
            assert_eq!(result.issues.len(), 1);
            let missing = result
                .warnings
                .iter()
                .filter(|v| matches!(v, ResolutionWarning::ModificationNotFound { .. }))
                .count();
            assert_eq!(missing, if mode == Mode::Try { 2 } else { 1 });
        }
        let mut original = db();
        assert!(
            run(
                Input::Chain(&pf),
                &mut original,
                Mode::Strict,
                None,
                &mut Budget::default()
            )
            .is_err()
        );
        assert!(original.is_empty());
    }
    #[test]
    fn shared_chain_work_cannot_reset_at_a_later_chain() {
        let chain = Peptidoform::parse("M[Formula:Cl101]PEPTIDE").unwrap();
        let mut one = Budget::default();
        run(
            Input::Chain(&chain),
            &mut db(),
            Mode::Strict,
            None,
            &mut one,
        )
        .unwrap();
        let used = MAX_PROFORMA_MASS_WORK - one.work;
        let ion = PeptidoformIon {
            chains: vec![chain.clone(), chain],
            ..Default::default()
        };
        let mut limited = Budget {
            work: used,
            ..Budget::default()
        };
        let mut original = db();
        assert!(
            run(
                Input::Ion(&ion),
                &mut original,
                Mode::Strict,
                None,
                &mut limited
            )
            .is_err()
        );
        assert!(original.is_empty());
    }
    #[test]
    fn early_empty_charge_and_chimeric_checks_do_not_copy_ignored_ast() {
        let mut none = Budget {
            work: 0,
            bytes: 0,
            items: 0,
        };
        assert_eq!(
            run(
                Input::Ion(&PeptidoformIon::default()),
                &mut db(),
                Mode::Try,
                None,
                &mut none
            )
            .unwrap()
            .value,
            Some(0.)
        );
        let chain = Peptidoform {
            name: Some("x".repeat(MAX_PROFORMA_MASS_TEXT_BYTES + 1)),
            ..Default::default()
        };
        let ion = PeptidoformIon {
            chains: vec![chain],
            is_chimeric: true,
            ..Default::default()
        };
        assert!(ion.try_mono_mass(&mut db()).unwrap().value.is_none());
        assert!(ion.try_mz(&mut db()).unwrap().value.is_none());
        assert!(ion.mass_calculation_issues(&mut db()).is_err());
    }
}
