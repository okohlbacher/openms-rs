// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::resolution::{Budget, Resolver, canonical_formula};
use super::*;
use crate::chemistry::{
    AASequence, EmpiricalFormula, MassTag, ModificationProvenance, ModificationRecord,
    ModificationsDB, SequenceModification, TermSpecificity, composition_formula,
    residue_composition,
};
use std::mem::size_of;

pub const MAX_PROFORMA_CONVERSION_WORK: usize = MAX_PROFORMA_RESOLUTION_WORK;
pub const MAX_PROFORMA_CONVERSION_BYTES: usize = MAX_PROFORMA_RESOLUTION_BYTES;
pub const MAX_PROFORMA_CONVERSION_ITEMS: usize = MAX_PROFORMA_RESOLUTION_ITEMS;
pub const MAX_PROFORMA_CONVERSION_TEXT_BYTES: usize = MAX_PROFORMA_RESOLUTION_TEXT_BYTES;

/// Owned source resolution and combination warnings, in emission order.
#[derive(Clone, Debug, PartialEq)]
pub enum ConversionWarning {
    Resolution(ResolutionWarning),
    FormulaMassDisagreement {
        formula: String,
        formula_mass: f64,
        summed_mass: f64,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConversionEvaluation<T> {
    pub value: T,
    pub warnings: Vec<ConversionWarning>,
}

impl Peptidoform {
    /// Source conversion diagnostics; completion can intern formula records.
    pub fn aa_sequence_conversion_issues(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<Vec<ConversionIssue>>> {
        conversion_issues(self, registry, &mut Budget::default())
    }
    pub fn is_representable_as_aa_sequence(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<bool>> {
        let result = self.aa_sequence_conversion_issues(registry)?;
        Ok(ConversionEvaluation {
            value: result.value.is_empty(),
            warnings: result.warnings,
        })
    }
    /// Convert with the source default FailOnLoss policy.
    pub fn to_aa_sequence(
        &self,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<AASequence>> {
        self.to_aa_sequence_with_policy(ConversionPolicy::FailOnLoss, registry)
    }
    /// Convert without changing the AST; any error also rolls back the registry.
    pub fn to_aa_sequence_with_policy(
        &self,
        policy: ConversionPolicy,
        registry: &mut ModificationsDB,
    ) -> Result<ConversionEvaluation<AASequence>> {
        convert(self, policy, registry, &mut Budget::default())
    }
    /// Convert immutable owned sequence annotations into source-shaped AST tags.
    /// Native anonymous MassTags obtain detached owned chemistry handles.
    pub fn from_aa_sequence(sequence: &AASequence) -> Result<Self> {
        from_sequence(sequence, &mut Budget::default())
    }
}

fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("ProForma conversion mass overflows"))
    }
}
fn vector<T>(length: usize, budget: &mut Budget) -> Result<Vec<T>> {
    budget.items(length)?;
    budget.allocate(length.saturating_mul(size_of::<T>()))?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(length)
        .map_err(|_| invalid("ProForma conversion allocation failed"))?;
    Ok(result)
}
fn copy_text(text: &str, budget: &mut Budget) -> Result<String> {
    budget.text(text)?;
    budget.allocate(text.len())?;
    Ok(text.to_owned())
}
fn finish<T>(
    value: T,
    resolver: Resolver<'_>,
    extra: Vec<ConversionWarning>,
) -> Result<(ConversionEvaluation<T>, Option<ModificationsDB>)> {
    // Precharge the combined vector before releasing the session's ledger.
    let total = resolver.warning_count().saturating_add(extra.len());
    let mut warnings = vector(total, resolver.budget)?;
    let (staged, resolution) = resolver.finish();
    warnings.extend(resolution.into_iter().map(ConversionWarning::Resolution));
    warnings.extend(extra);
    Ok((ConversionEvaluation { value, warnings }, staged))
}
fn conversion_issues(
    input: &Peptidoform,
    registry: &mut ModificationsDB,
    budget: &mut Budget,
) -> Result<ConversionEvaluation<Vec<ConversionIssue>>> {
    let mut resolved = super::mass::clone_chain(input, budget)?;
    let mut resolver = Resolver::new(registry, budget);
    resolver.resolve_chain(&mut resolved)?;
    let issues = collect_issues(&resolved, resolver.budget)?;
    let (result, staged) = finish(issues, resolver, Vec::new())?;
    if let Some(next) = staged {
        *registry = next;
    }
    Ok(result)
}

fn chemistry(tag: &ModificationTag) -> bool {
    !matches!(
        tag,
        ModificationTag::InfoTag(_) | ModificationTag::PositionConstraint(_)
    )
}
fn chemistry_count(modification: &Modification, budget: &mut Budget) -> Result<usize> {
    budget.items(modification.alternatives.len())?;
    Ok(modification
        .alternatives
        .iter()
        .filter(|(tag, _)| chemistry(tag))
        .count())
}
fn issue(
    out: &mut Vec<ConversionIssue>,
    kind: ConversionIssueType,
    position: Option<usize>,
    description: impl FnOnce() -> String,
    budget: &mut Budget,
) -> Result<()> {
    // Every source message below has bounded literals plus at most a usize.
    budget.consume(1024)?;
    budget.allocate(1024)?;
    budget.reserve(out)?;
    out.push(ConversionIssue {
        issue_type: kind,
        description: description(),
        position,
    });
    Ok(())
}
fn collect_issues(pf: &Peptidoform, budget: &mut Budget) -> Result<Vec<ConversionIssue>> {
    let mut out = Vec::new();
    for (present, kind, text) in [
        (
            !pf.unlocalised_mods.is_empty(),
            ConversionIssueType::UnlocalisedMod,
            "Peptidoform contains unlocalised modifications",
        ),
        (
            !pf.labile_mods.is_empty(),
            ConversionIssueType::LabileMod,
            "Peptidoform contains labile modifications",
        ),
    ] {
        if present {
            issue(&mut out, kind, None, || text.into(), budget)?;
        }
    }
    budget.items(pf.global_mods.len())?;
    if pf
        .global_mods
        .iter()
        .any(|entry| matches!(entry, GlobalModEntry::GlobalModification(_)))
    {
        issue(
            &mut out,
            ConversionIssueType::GlobalMod,
            None,
            || "Peptidoform contains global modifications".into(),
            budget,
        )?;
    }
    budget.items(pf.sequence.len())?;
    for (position, section) in pf.sequence.iter().enumerate() {
        match section {
            SequenceSection::AmbiguousRegion(_) => issue(
                &mut out,
                ConversionIssueType::AmbiguousRegion,
                Some(position),
                || format!("Peptidoform contains ambiguous region at position {position}"),
                budget,
            )?,
            SequenceSection::ModifiedRange(_) => issue(
                &mut out,
                ConversionIssueType::ModifiedRange,
                Some(position),
                || format!("Peptidoform contains modified range at position {position}"),
                budget,
            )?,
            SequenceSection::Element(element) => {
                budget.items(element.modifications.len())?;
                let mut brackets = 0;
                for modification in &element.modifications {
                    brackets += usize::from(chemistry_count(modification, budget)? != 0);
                }
                if brackets > 1 {
                    issue(
                        &mut out,
                        ConversionIssueType::UnsupportedFeature,
                        Some(position),
                        || {
                            format!(
                                "Residue at position {position} carries multiple modification brackets; an AASequence residue holds only one, so they are combined into one anonymous modification and the individual identities are lost"
                            )
                        },
                        budget,
                    )?;
                }
                for modification in &element.modifications {
                    let count = chemistry_count(modification, budget)?;
                    if modification.resolved_mod.is_none() && count != 0 {
                        issue(
                            &mut out,
                            ConversionIssueType::UnresolvedMod,
                            Some(position),
                            || format!("Modification at position {position} could not be resolved"),
                            budget,
                        )?;
                    }
                    if count > 1 {
                        issue(
                            &mut out,
                            ConversionIssueType::AlternativeMods,
                            Some(position),
                            || {
                                format!(
                                    "Modification at position {position} has multiple alternatives"
                                )
                            },
                            budget,
                        )?;
                    }
                    budget.items(modification.alternatives.len())?;
                    if modification.alternatives.iter().any(|(_, label)| {
                        label
                            .as_ref()
                            .is_some_and(|label| label.label_type == LabelType::Crosslink)
                    }) {
                        issue(
                            &mut out,
                            ConversionIssueType::CrossLink,
                            Some(position),
                            || {
                                format!(
                                    "Modification at position {position} is part of a cross-link"
                                )
                            },
                            budget,
                        )?;
                    }
                }
            }
        }
    }
    for (groups, label) in [(&pf.n_term_mods, "N"), (&pf.c_term_mods, "C")] {
        budget.items(groups.len())?;
        let mut brackets = 0;
        for modification in groups {
            brackets += usize::from(chemistry_count(modification, budget)? != 0);
        }
        if brackets > 1 {
            issue(
                &mut out,
                ConversionIssueType::UnsupportedFeature,
                None,
                || {
                    format!(
                        "{label}-terminal modifications: an AASequence holds only one, so all but the first are lost"
                    )
                },
                budget,
            )?;
        }
        for modification in groups {
            let count = chemistry_count(modification, budget)?;
            if modification.resolved_mod.is_none() && count != 0 {
                issue(
                    &mut out,
                    ConversionIssueType::UnresolvedMod,
                    None,
                    || format!("{label}-terminal modification could not be resolved"),
                    budget,
                )?;
            }
            if count > 1 {
                issue(
                    &mut out,
                    ConversionIssueType::AlternativeMods,
                    None,
                    || format!("{label}-terminal modification has multiple alternatives"),
                    budget,
                )?;
            }
        }
    }
    Ok(out)
}

fn convert(
    input: &Peptidoform,
    policy: ConversionPolicy,
    registry: &mut ModificationsDB,
    budget: &mut Budget,
) -> Result<ConversionEvaluation<AASequence>> {
    let mut pf = super::mass::clone_chain(input, budget)?;
    let mut resolver = Resolver::new(registry, budget);
    resolver.resolve_chain(&mut pf)?;
    let issues = collect_issues(&pf, resolver.budget)?;
    if policy == ConversionPolicy::FailOnLoss && !issues.is_empty() {
        let length = issues.iter().fold(45usize, |n, issue| {
            n.saturating_add(issue.description.len()).saturating_add(2)
        });
        if length > MAX_PROFORMA_CONVERSION_TEXT_BYTES {
            return Err(invalid(
                "ProForma conversion diagnostic text limit exceeded",
            ));
        }
        resolver.budget.consume(length)?;
        resolver.budget.allocate(length)?;
        let mut message = String::with_capacity(length);
        message.push_str("Cannot convert Peptidoform to AASequence: ");
        for issue in issues {
            message.push_str(&issue.description);
            message.push_str("; ");
        }
        return Err(invalid(message));
    }
    resolver.budget.items(pf.sequence.len())?;
    let mut chars = 0usize;
    for section in &pf.sequence {
        chars = chars.saturating_add(match section {
            SequenceSection::Element(_) => 1,
            SequenceSection::AmbiguousRegion(region) => usize::from(!region.elements.is_empty()),
            SequenceSection::ModifiedRange(range) => range.elements.len(),
        });
    }
    resolver.budget.items(chars)?;
    resolver.budget.allocate(chars.saturating_mul(4))?;
    let mut text = String::new();
    text.try_reserve_exact(chars.saturating_mul(4))
        .map_err(|_| invalid("ProForma conversion allocation failed"))?;
    for section in &pf.sequence {
        match section {
            SequenceSection::Element(element) => text.push(element.amino_acid),
            SequenceSection::AmbiguousRegion(region) => {
                if let Some(element) = region.elements.first() {
                    text.push(element.amino_acid);
                }
            }
            SequenceSection::ModifiedRange(range) => {
                for element in &range.elements {
                    text.push(element.amino_acid);
                }
            }
        }
    }
    resolver.budget.text(&text)?;
    let (mut work, mut bytes) = (resolver.budget.work, resolver.budget.bytes);
    let base = AASequence::parse_with_budget(&text, resolver.database(), &mut work, &mut bytes);
    resolver.budget.work = work;
    resolver.budget.bytes = bytes;
    let base = base?;
    let mut assignments =
        vector::<(usize, Arc<ResidueModification>)>(pf.sequence.len(), resolver.budget)?;
    let mut warnings = Vec::new();
    let mut position = 0usize;
    for section in &pf.sequence {
        match section {
            SequenceSection::Element(element) => {
                let mut resolved = vector(element.modifications.len(), resolver.budget)?;
                for modification in &element.modifications {
                    if let Some(record) = &modification.resolved_mod {
                        resolved.push(Arc::clone(record));
                    } else if policy == ConversionPolicy::FailOnLoss {
                        resolver.budget.consume(128)?;
                        resolver.budget.allocate(128)?;
                        return Err(invalid(format!(
                            "Unresolved modification at position {position}"
                        )));
                    }
                }
                let record = match resolved.len() {
                    0 => None,
                    1 => resolved.pop(),
                    _ => combine(&resolved, element.amino_acid, &mut resolver, &mut warnings)?,
                };
                if let Some(record) = record {
                    // Source checks this index at each actual attachment, before
                    // later sections can perform additional work or interning.
                    if position >= base.len() {
                        return Err(invalid("resolved modification index out of bounds"));
                    }
                    assignments.push((position, record));
                }
                position = position
                    .checked_add(1)
                    .ok_or_else(|| invalid("conversion position overflows"))?;
            }
            SequenceSection::AmbiguousRegion(_) => {
                position = position
                    .checked_add(1)
                    .ok_or_else(|| invalid("conversion position overflows"))?;
            }
            SequenceSection::ModifiedRange(range) => {
                position = position
                    .checked_add(range.elements.len())
                    .ok_or_else(|| invalid("conversion position overflows"))?;
            }
        }
    }
    resolver
        .budget
        .items(pf.n_term_mods.len().saturating_add(pf.c_term_mods.len()))?;
    let n = pf
        .n_term_mods
        .iter()
        .find_map(|m| m.resolved_mod.as_ref())
        .cloned();
    let c = pf
        .c_term_mods
        .iter()
        .find_map(|m| m.resolved_mod.as_ref())
        .cloned();
    charge_attachment(
        &base,
        &assignments,
        n.as_deref(),
        c.as_deref(),
        resolver.budget,
    )?;
    let result = base.with_resolved_modifications(&assignments, n, c)?;
    let (result, staged) = finish(result, resolver, warnings)?;
    if let Some(next) = staged {
        *registry = next;
    }
    Ok(result)
}

fn base_formula(residue: char) -> Option<EmpiricalFormula> {
    if !residue.is_ascii() {
        return None;
    }
    residue_composition(residue as u8)
        .map(|c| {
            composition_formula(c)
                .checked_add(&water())
                .expect("small source residue formula")
        })
        .or_else(|| matches!(residue, 'B' | 'Z' | 'X').then(EmpiricalFormula::default))
}
fn water() -> EmpiricalFormula {
    composition_formula([0, 2, 0, 1, 0, 0])
}
fn combine(
    records: &[Arc<ResidueModification>],
    residue: char,
    resolver: &mut Resolver<'_>,
    warnings: &mut Vec<ConversionWarning>,
) -> Result<Option<Arc<ResidueModification>>> {
    if residue == '\0' {
        return Ok(records.last().cloned());
    }
    resolver.budget.items(records.len())?;
    let mut mass = 0.;
    let mut formula = EmpiricalFormula::default();
    let mut all_formulas = true;
    for record in records {
        resolver.budget.record(record)?;
        mass = finite(mass + record.diff_mono_mass())?;
        if record.diff_formula().is_empty() {
            all_formulas = false;
        } else {
            let atoms = formula
                .stored_atom_types()
                .saturating_add(record.diff_formula().stored_atom_types())
                .saturating_add(1);
            resolver.budget.consume(atoms.saturating_mul(256))?;
            resolver.budget.allocate(atoms.saturating_mul(1024))?;
            formula = formula.checked_add(record.diff_formula())?;
        }
    }
    if mass.abs() <= 1e-6 {
        return Ok(None);
    }
    resolver.budget.consume(
        formula
            .stored_atom_types()
            .saturating_add(1)
            .saturating_mul(128),
    )?;
    let formula_mass = finite(formula.mono_mass())?;
    if all_formulas && !formula.is_empty() && (formula_mass - mass).abs() <= 1e-3 {
        let tag = FormulaTag {
            formula_string: canonical_formula(&formula, resolver.budget)?,
            charge: None,
        };
        return Ok(resolver
            .formula(&tag, Some(residue), Some(TermSpecificity::Anywhere))?
            .or_else(|| records.last().cloned()));
    }
    if all_formulas {
        let text = canonical_formula(&formula, resolver.budget)?;
        resolver.budget.reserve(warnings)?;
        warnings.push(ConversionWarning::FormulaMassDisagreement {
            formula: text,
            formula_mass,
            summed_mass: mass,
        });
    }
    let Some(base) = base_formula(residue) else {
        return Ok(records.last().cloned());
    };
    resolver.budget.consume(4096)?;
    resolver.budget.allocate(4096)?;
    let number = signed_source_mass(mass);
    let full_id = format!("{residue}[{number}]");
    if let Some(existing) = resolver.existing_full_id(&full_id)? {
        return Ok(Some(existing));
    }
    let record = ResidueModification::from_record(ModificationRecord {
        full_id,
        full_name: format!("[{number}]"),
        origin: Some(residue),
        diff_mono_mass: mass,
        mono_mass: finite(mass + base.mono_mass())?,
        average_mass: finite(mass + base.average_mass())?,
        provenance: ModificationProvenance::MassOnly,
        ..ModificationRecord::default()
    })?;
    resolver.insert_record(record)
}
fn signed_source_mass(mass: f64) -> String {
    format!(
        "{}{}",
        if mass < 0. { "-" } else { "+" },
        crate::param::value::format_float(mass.abs(), true)
    )
}

fn charge_attachment(
    base: &AASequence,
    assignments: &[(usize, Arc<ResidueModification>)],
    n: Option<&ResidueModification>,
    c: Option<&ResidueModification>,
    budget: &mut Budget,
) -> Result<()> {
    budget.items(base.len().saturating_add(assignments.len()))?;
    let payload = base.generation_payload_bytes()?;
    budget.consume(payload)?;
    budget.allocate(payload)?;
    let mut atoms = 6usize;
    for record in assignments
        .iter()
        .map(|(_, m)| m.as_ref())
        .chain(n)
        .chain(c)
    {
        budget.record(record)?;
        atoms = atoms
            .saturating_add(record.diff_formula().stored_atom_types())
            .saturating_add(
                record
                    .absolute_formula()
                    .map_or(0, EmpiricalFormula::stored_atom_types),
            );
    }
    // Existing AASequence rebuild clones cumulative formula maps; conservatively
    // precharge every residue against the union-size upper bound before calling it.
    let nodes = base
        .len()
        .saturating_add(3)
        .saturating_mul(atoms.saturating_add(12));
    budget.consume(nodes.saturating_mul(256))?;
    budget.allocate(nodes.saturating_mul(1024))?;
    Ok(())
}

fn from_sequence(sequence: &AASequence, budget: &mut Budget) -> Result<Peptidoform> {
    let mut pf = Peptidoform {
        sequence: vector(sequence.len(), budget)?,
        ..Peptidoform::default()
    };
    for (index, residue) in sequence.as_str().chars().enumerate() {
        let modifications = match sequence.residue_modification(index)? {
            Some(modification) => {
                let mut group = vector(1, budget)?;
                group.push(from_modification(modification, budget)?);
                group
            }
            None => Vec::new(),
        };
        pf.sequence.push(SequenceSection::Element(SequenceElement {
            amino_acid: residue,
            modifications,
        }));
    }
    for (modification, group) in [
        (sequence.n_terminal_modification(), &mut pf.n_term_mods),
        (sequence.c_terminal_modification(), &mut pf.c_term_mods),
    ] {
        if let Some(modification) = modification {
            *group = vector(1, budget)?;
            group.push(from_modification(modification, budget)?);
        }
    }
    Ok(pf)
}
fn from_modification(
    modification: &SequenceModification,
    budget: &mut Budget,
) -> Result<Modification> {
    let record = match modification {
        SequenceModification::Known(record) => {
            budget.record(record)?;
            Arc::clone(record)
        }
        SequenceModification::MassTag(tag) => mass_tag_record(tag, budget)?,
    };
    let mut alternatives = vector(2, budget)?;
    if let Some(id) = record.record_id() {
        budget.consume(64)?;
        budget.allocate(64)?;
        alternatives.push((
            ModificationTag::CvAccession(CvAccession {
                database: CvDatabase::Unimod,
                accession: id.to_string(),
            }),
            None,
        ));
    } else if record.name().is_empty() && !record.diff_formula().is_empty() {
        alternatives.push((formula_tag(record.diff_formula(), budget)?, None));
    } else if !record.name().is_empty() && record.provenance() == ModificationProvenance::Defined {
        let tag = if !record.diff_formula().is_empty() {
            formula_tag(record.diff_formula(), budget)?
        } else {
            mass_tag(record.diff_mono_mass(), budget)?
        };
        alternatives.push((tag, None));
        alternatives.push((
            ModificationTag::InfoTag(InfoTag {
                text: copy_text(record.name(), budget)?,
            }),
            None,
        ));
    } else if !record.name().is_empty() {
        alternatives.push((
            ModificationTag::NamedMod(NamedMod {
                name: copy_text(record.name(), budget)?,
                cv_hint: None,
            }),
            None,
        ));
    } else {
        alternatives.push((mass_tag(record.diff_mono_mass(), budget)?, None));
    }
    Ok(Modification {
        alternatives,
        resolved_mod: Some(record),
    })
}
fn formula_tag(formula: &EmpiricalFormula, budget: &mut Budget) -> Result<ModificationTag> {
    Ok(ModificationTag::FormulaTag(FormulaTag {
        formula_string: canonical_formula(formula, budget)?,
        charge: None,
    }))
}
fn mass_tag(mass: f64, budget: &mut Budget) -> Result<ModificationTag> {
    finite(mass)?;
    budget.consume(2048)?;
    budget.allocate(2048)?;
    // Integral doubles require exact fixed zero precision: shortest general
    // Display may round low integer digits (e.g. 100000000000000016384).
    // For fractional doubles Display supplies shortest fixed decimal. The
    // extracted source formatter fixture checks both cases and subnormals.
    let magnitude = if mass.fract() == 0. {
        format!("{:.0}", mass.abs())
    } else {
        mass.abs().to_string()
    };
    let original_text = format!("{}{}", if mass < 0. { "-" } else { "+" }, magnitude);
    Ok(ModificationTag::MassDelta(MassDelta {
        mass,
        source: MassDeltaSource::None,
        original_text,
    }))
}
fn mass_tag_record(tag: &MassTag, budget: &mut Budget) -> Result<Arc<ResidueModification>> {
    budget.text(tag.full_id())?;
    budget.text(tag.input())?;
    budget.consume(
        tag.full_id()
            .len()
            .saturating_add(tag.input().len())
            .saturating_add(4096),
    )?;
    budget.allocate(
        tag.full_id()
            .len()
            .saturating_add(tag.input().len())
            .saturating_add(4096),
    )?;
    let term = tag.term_specificity();
    let base = if term == TermSpecificity::Anywhere {
        base_formula(
            tag.origin()
                .ok_or_else(|| invalid("mass tag has no origin"))?,
        )
        .ok_or_else(|| invalid("mass tag has unknown origin"))?
    } else if matches!(term, TermSpecificity::NTerm | TermSpecificity::ProteinNTerm) {
        composition_formula([0, 1, 0, 0, 0, 0])
    } else {
        composition_formula([0, 1, 0, 1, 0, 0])
    };
    let delta = finite(match tag.delta_mono_mass() {
        Some(delta) => delta,
        // Only native absolute B/Z/X tags lack a delta; source full mass is zero.
        None => {
            tag.residue_mono_mass()
                .ok_or_else(|| invalid("mass tag has no target mass"))?
                + water().mono_mass()
        }
    })?;
    let absolute = finite(delta + base.mono_mass())?;
    let average = if term == TermSpecificity::Anywhere {
        finite(delta + base.average_mass())?
    } else {
        0.
    };
    Ok(Arc::new(ResidueModification::from_record(
        ModificationRecord {
            full_id: tag.full_id().to_owned(),
            full_name: format!("[{}]", tag.input()),
            origin: tag.origin(),
            term_specificity: term,
            diff_mono_mass: delta,
            mono_mass: absolute,
            average_mass: average,
            provenance: ModificationProvenance::MassOnly,
            ..ModificationRecord::default()
        },
    )?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn late_failures_after_interning_roll_back_every_publication() {
        let pf = Peptidoform::parse("M[Formula:Cl101]A").unwrap();
        let snapshot = pf.clone();
        let mut complete = Budget::default();
        let mut db = ModificationsDB::default();
        convert(&pf, ConversionPolicy::FailOnLoss, &mut db, &mut complete).unwrap();
        assert_eq!(db.len(), 1);
        let used = [
            MAX_PROFORMA_CONVERSION_WORK - complete.work,
            MAX_PROFORMA_CONVERSION_BYTES - complete.bytes,
            MAX_PROFORMA_CONVERSION_ITEMS - complete.items,
        ];
        for (field, limit) in used.into_iter().enumerate() {
            let mut budget = Budget::default();
            match field {
                0 => budget.work = limit - 1,
                1 => budget.bytes = limit - 1,
                _ => budget.items = limit - 1,
            }
            let mut db = ModificationsDB::default();
            assert!(convert(&pf, ConversionPolicy::FailOnLoss, &mut db, &mut budget).is_err());
            assert!(db.is_empty());
            assert_eq!(pf, snapshot);
        }
    }

    #[test]
    fn shared_copy_and_attachment_costs_do_not_reset_for_later_residues() {
        let one = Peptidoform::parse("M[Formula:Cl101]").unwrap();
        let two = Peptidoform::parse("M[Formula:Cl101]M[Formula:Cl101]").unwrap();
        let mut used = Budget::default();
        convert(
            &one,
            ConversionPolicy::BestEffort,
            &mut ModificationsDB::default(),
            &mut used,
        )
        .unwrap();
        let mut budget = Budget {
            work: MAX_PROFORMA_CONVERSION_WORK - used.work,
            ..Budget::default()
        };
        let mut db = ModificationsDB::default();
        assert!(convert(&two, ConversionPolicy::BestEffort, &mut db, &mut budget).is_err());
        assert!(db.is_empty());
        let mut zero = Budget {
            bytes: 0,
            ..Budget::default()
        };
        assert!(from_sequence(&AASequence::parse("M").unwrap(), &mut zero).is_err());
    }

    #[test]
    fn reverse_fixed_decimal_extremes_are_checked_and_roundtrip() {
        for mass in [
            f64::MAX,
            -f64::MAX,
            f64::from_bits(1),
            -f64::from_bits(1),
            0.,
            -0.,
        ] {
            let tag = mass_tag(mass, &mut Budget::default()).unwrap();
            let ModificationTag::MassDelta(delta) = tag else {
                panic!()
            };
            assert!(!delta.original_text.contains(['e', 'E']));
            assert!(delta.original_text.len() <= 400);
            let parsed: f64 = delta.original_text.parse().unwrap();
            assert_eq!(parsed, mass);
            if mass != 0. {
                assert_eq!(parsed.to_bits(), mass.to_bits());
            } else {
                assert_eq!(delta.original_text, "+0");
            }
        }
        assert!(mass_tag(f64::INFINITY, &mut Budget::default()).is_err());
        let mut budget = Budget {
            bytes: 2047,
            ..Budget::default()
        };
        assert!(mass_tag(1., &mut budget).is_err());
    }

    #[test]
    fn sparse_formula_and_record_work_is_precharged_before_rebuild() {
        let base = AASequence::parse("M").unwrap();
        let record = Arc::new(
            ResidueModification::from_record(ModificationRecord {
                name: "custom".into(),
                origin: Some('M'),
                diff_formula: EmpiricalFormula::parse("Cl").unwrap(),
                ..ModificationRecord::default()
            })
            .unwrap(),
        );
        let assignments = [(0, record)];
        let mut budget = Budget {
            bytes: 1024,
            ..Budget::default()
        };
        assert!(charge_attachment(&base, &assignments, None, None, &mut budget).is_err());
        assert!(!base.is_modified());
    }
}
