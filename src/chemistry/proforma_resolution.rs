// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::*;
use crate::chemistry::{
    EmpiricalFormula, ModificationProvenance, ModificationRecord, ModificationsDB, TermSpecificity,
    composition_formula, residue_composition,
};
use std::mem::size_of;

pub const MAX_PROFORMA_RESOLUTION_WORK: usize = 50_000_000;
pub const MAX_PROFORMA_RESOLUTION_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_PROFORMA_RESOLUTION_ITEMS: usize = 1_000_000;
pub const MAX_PROFORMA_RESOLUTION_TEXT_BYTES: usize = 4 * 1024 * 1024;

/// Owned equivalents of warnings reachable through source modification resolution.
/// Lookup ambiguity selects first native provider order instead of pointer order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ResolutionWarning {
    ModificationNotFound {
        name: String,
    },
    AmbiguousAccession {
        name: String,
        residue: Option<char>,
        term: Option<TermSpecificity>,
        selected: Arc<ResidueModification>,
    },
    DefinitionDisagreement {
        definition: Arc<ResidueModification>,
        inline: Arc<ResidueModification>,
    },
}

impl Peptidoform {
    /// Resolve the source-supported annotation sites against a caller-owned DB.
    /// New formula definitions and handle replacements are published together;
    /// checked failure leaves both objects unchanged. Empty alternatives retain
    /// their old handle. See `docs/PROFORMA_RESOLUTION_SUPPORT.md` for traversal,
    /// lookup ordering and the source's intentionally unresolved tag forms.
    pub fn resolve_modifications(
        &mut self,
        registry: &mut ModificationsDB,
    ) -> Result<Vec<ResolutionWarning>> {
        resolve(self, registry, &mut Budget::default())
    }
}

pub(super) struct Budget {
    pub(super) work: usize,
    pub(super) bytes: usize,
    pub(super) items: usize,
}
impl Default for Budget {
    fn default() -> Self {
        Self {
            work: MAX_PROFORMA_RESOLUTION_WORK,
            bytes: MAX_PROFORMA_RESOLUTION_BYTES,
            items: MAX_PROFORMA_RESOLUTION_ITEMS,
        }
    }
}
impl Budget {
    pub(super) fn consume(&mut self, amount: usize) -> Result<()> {
        self.work = self.work.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    pub(super) fn allocate(&mut self, amount: usize) -> Result<()> {
        self.bytes = self.bytes.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    pub(super) fn items(&mut self, amount: usize) -> Result<()> {
        self.items = self.items.checked_sub(amount).ok_or_else(limit)?;
        self.consume(amount)
    }
    pub(super) fn text(&mut self, text: &str) -> Result<()> {
        if text.len() > MAX_PROFORMA_RESOLUTION_TEXT_BYTES {
            return Err(limit());
        }
        self.consume(text.len().saturating_add(1))
    }
    pub(super) fn reserve<T>(&mut self, values: &mut Vec<T>) -> Result<()> {
        if values.len() == values.capacity() {
            let additional = values.len().max(1);
            self.consume(values.len().saturating_add(additional))?;
            self.allocate(
                values
                    .len()
                    .saturating_add(additional)
                    .saturating_mul(size_of::<T>()),
            )?;
            values.try_reserve_exact(additional).map_err(|_| limit())?;
        }
        Ok(())
    }
    pub(super) fn record(&mut self, record: &ResidueModification) -> Result<()> {
        // Precharge collection walks before payload_bytes visits arbitrary
        // caller-defined synonyms/losses. Formula maps themselves need only len.
        self.consume(
            record
                .synonyms()
                .len()
                .saturating_add(record.neutral_losses().len())
                .saturating_add(16),
        )?;
        let mut atoms = record.diff_formula().stored_atom_types().saturating_add(
            record
                .absolute_formula()
                .map_or(0, EmpiricalFormula::stored_atom_types),
        );
        for loss in record.neutral_losses() {
            atoms = atoms.saturating_add(loss.formula().stored_atom_types());
        }
        self.consume(atoms.saturating_add(record.payload_bytes()?))
    }
}
fn limit() -> Error {
    invalid("ProForma resolution resource limit exceeded")
}

type Patch<'a> = (&'a mut Modification, Option<Arc<ResidueModification>>);

fn resolve(
    pf: &mut Peptidoform,
    registry: &mut ModificationsDB,
    budget: &mut Budget,
) -> Result<Vec<ResolutionWarning>> {
    let mut resolver = Resolver::new(registry, budget);
    resolver.resolve_chain(pf)?;
    let (staged, warnings) = resolver.finish();
    if let Some(next) = staged {
        *registry = next;
    }
    Ok(warnings)
}

pub(super) struct Resolver<'a> {
    original: &'a ModificationsDB,
    staged: Option<ModificationsDB>,
    pub(super) budget: &'a mut Budget,
    warnings: Vec<ResolutionWarning>,
}
impl<'a> Resolver<'a> {
    pub(super) fn new(original: &'a ModificationsDB, budget: &'a mut Budget) -> Self {
        Self {
            original,
            staged: None,
            budget,
            warnings: Vec::new(),
        }
    }
    pub(super) fn warning_count(&self) -> usize {
        self.warnings.len()
    }
    pub(super) fn finish(self) -> (Option<ModificationsDB>, Vec<ResolutionWarning>) {
        (self.staged, self.warnings)
    }
    /// Stage one chain's complete handle patch before assigning any handle.
    /// The registry stays in this private session until its caller publishes it.
    pub(super) fn resolve_chain(&mut self, pf: &mut Peptidoform) -> Result<()> {
        let mut patches = Vec::new();
        self.budget.items(pf.sequence.len())?;
        for section in &mut pf.sequence {
            match section {
                SequenceSection::Element(element) => {
                    self.group(
                        &mut element.modifications,
                        Some(element.amino_acid),
                        None,
                        &mut patches,
                    )?;
                }
                SequenceSection::AmbiguousRegion(region) => {
                    self.budget.items(region.elements.len())?;
                    for element in &mut region.elements {
                        self.group(
                            &mut element.modifications,
                            Some(element.amino_acid),
                            None,
                            &mut patches,
                        )?;
                    }
                }
                SequenceSection::ModifiedRange(range) => {
                    // Literal source omission: modifications on range.elements
                    // are not traversed by resolveModifications.
                    self.group(&mut range.modifications, None, None, &mut patches)?;
                }
            }
        }
        self.group(
            &mut pf.n_term_mods,
            None,
            Some(TermSpecificity::NTerm),
            &mut patches,
        )?;
        self.group(
            &mut pf.c_term_mods,
            None,
            Some(TermSpecificity::CTerm),
            &mut patches,
        )?;
        self.budget.items(pf.unlocalised_mods.len())?;
        for group in &mut pf.unlocalised_mods {
            self.group(&mut group.modifications, None, None, &mut patches)?;
        }
        self.budget.items(pf.labile_mods.len())?;
        for labile in &mut pf.labile_mods {
            self.group(
                std::slice::from_mut(&mut labile.modification),
                None,
                None,
                &mut patches,
            )?;
        }
        self.budget.items(pf.global_mods.len())?;
        for entry in &mut pf.global_mods {
            if let GlobalModEntry::GlobalModification(global) = entry {
                self.group(
                    std::slice::from_mut(&mut global.modification),
                    None,
                    None,
                    &mut patches,
                )?;
            }
        }
        for (target, value) in patches {
            target.resolved_mod = value;
        }
        Ok(())
    }
}
impl Resolver<'_> {
    pub(super) fn database(&self) -> &ModificationsDB {
        self.staged.as_ref().unwrap_or(self.original)
    }
    fn group<'a>(
        &mut self,
        modifications: &'a mut [Modification],
        residue: Option<char>,
        term: Option<TermSpecificity>,
        patches: &mut Vec<Patch<'a>>,
    ) -> Result<()> {
        let residue = residue.filter(|&value| value != '\0');
        self.budget.items(modifications.len())?;
        for modification in modifications {
            if modification.alternatives.is_empty() {
                continue;
            }
            self.budget.items(modification.alternatives.len())?;
            self.budget
                .consume(modification.alternatives.len().saturating_mul(2))?;
            let value = self.modification(modification, residue, term)?;
            if let Some(old) = &modification.resolved_mod {
                self.budget.record(old)?;
            }
            self.budget.reserve(patches)?;
            patches.push((modification, value));
        }
        Ok(())
    }
    fn warn(&mut self, warning: ResolutionWarning) -> Result<()> {
        self.budget.items(1)?;
        self.budget.reserve(&mut self.warnings)?;
        self.warnings.push(warning);
        Ok(())
    }
    fn lookup_charge(&mut self, name: &str) -> Result<()> {
        self.budget.text(name)?;
        // A usize-sized BTree needs fewer than1024 comparisons. Query length
        // bounds each comparison even for very long unrelated stored keys.
        self.budget.consume(
            name.len()
                .saturating_add(1)
                .saturating_mul(1024)
                .saturating_add(self.database().len().saturating_mul(4)),
        )
    }
    fn lookup(
        &mut self,
        name: &str,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<(Option<Arc<ResidueModification>>, usize)> {
        self.lookup_charge(name)?;
        let mut result = select_exact(self.database(), name, residue, term, false);
        let mut normalized = None;
        if !result.0
            && name.len() > 6
            && name
                .get(..6)
                .is_some_and(|prefix| prefix.eq_ignore_ascii_case("unimod"))
        {
            self.budget.allocate(name.len())?;
            self.budget
                .consume(name.len().saturating_add(1).saturating_mul(1024))?;
            let value = format!("UniMod{}", &name[6..]);
            result = select_exact(self.database(), &value, residue, term, false);
            normalized = Some(value);
        }
        if !result.0 {
            self.budget.allocate(name.len())?;
            self.budget.consume(name.len())?;
            self.warn(ResolutionWarning::ModificationNotFound {
                name: normalized.unwrap_or_else(|| name.into()),
            })?;
        }
        Ok((result.1, result.2))
    }
    fn defined(&mut self, name: &str) -> Result<bool> {
        self.lookup_charge(name)?;
        // Source gate uses exact aliases, not normalized UniMod spelling.
        Ok(select_exact(self.database(), name, None, None, true)
            .1
            .is_some())
    }
    fn modification(
        &mut self,
        modification: &Modification,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Option<Arc<ResidueModification>>> {
        let mut by_name = None;
        for (tag, _) in &modification.alternatives {
            if let ModificationTag::InfoTag(info) = tag {
                if !info.text.is_empty() && self.defined(&info.text)? {
                    by_name = self.lookup(&info.text, residue, term)?.0;
                    if by_name.is_some() {
                        break;
                    }
                }
            }
        }
        let chemistry = modification.alternatives.iter().find(|(tag, _)| {
            !matches!(
                tag,
                ModificationTag::InfoTag(_) | ModificationTag::PositionConstraint(_)
            )
        });
        let by_tag = match chemistry {
            Some((tag, _)) => self.tag(tag, residue, term)?,
            None => None,
        };
        if let Some(definition) = by_name {
            if let Some(inline) = by_tag {
                self.budget.consume(
                    definition
                        .diff_formula()
                        .stored_atom_types()
                        .saturating_add(inline.diff_formula().stored_atom_types())
                        .saturating_add(1),
                )?;
                let same =
                    if !definition.diff_formula().is_empty() && !inline.diff_formula().is_empty() {
                        definition.diff_formula() == inline.diff_formula()
                    } else {
                        (definition.diff_mono_mass() - inline.diff_mono_mass()).abs() <= 1e-6
                    };
                if !same {
                    self.warn(ResolutionWarning::DefinitionDisagreement {
                        definition: Arc::clone(&definition),
                        inline,
                    })?;
                }
            }
            return Ok(Some(definition));
        }
        if chemistry.is_some() {
            return Ok(by_tag);
        }
        self.tag(&modification.alternatives[0].0, residue, term)
    }
    fn tag(
        &mut self,
        tag: &ModificationTag,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Option<Arc<ResidueModification>>> {
        let residue = residue.filter(|&value| value != '\0');
        match tag {
            ModificationTag::CvAccession(accession) => {
                self.budget.text(&accession.accession)?;
                self.budget
                    .allocate(accession.accession.len().saturating_add(8))?;
                let prefix = match accession.database {
                    CvDatabase::Unimod => "UNIMOD",
                    CvDatabase::Mod => "MOD",
                    CvDatabase::Resid => "RESID",
                    CvDatabase::Xlmod => "XLMOD",
                    CvDatabase::Gno => "GNO",
                };
                let name = format!("{prefix}:{}", accession.accession);
                let mut selected = (None, 0);
                if residue.is_some() && term.is_none() {
                    selected = self.lookup(&name, residue, Some(TermSpecificity::Anywhere))?;
                }
                if selected.0.is_none() {
                    selected = self.lookup(&name, residue, term)?;
                }
                if selected.1 > 1 {
                    self.warn(ResolutionWarning::AmbiguousAccession {
                        name,
                        residue,
                        term,
                        selected: Arc::clone(selected.0.as_ref().expect("multiple matches")),
                    })?;
                }
                Ok(selected.0)
            }
            ModificationTag::NamedMod(named) => Ok(self.lookup(&named.name, residue, term)?.0),
            ModificationTag::MassDelta(delta) => {
                if !delta.mass.is_finite() {
                    return Err(invalid("ProForma resolution mass must be finite"));
                }
                self.budget
                    .consume(self.database().len().saturating_mul(8).saturating_add(1))?;
                let mut best = None;
                let mut error = 0.01;
                for record in self.database().entries() {
                    let difference = (record.diff_mono_mass() - delta.mass).abs();
                    if difference < error
                        && (record.diff_mono_mass() != 0. || delta.mass == 0.)
                        && source_matches(record, residue, term)
                    {
                        best = Some(Arc::clone(record));
                        error = difference;
                    }
                }
                Ok(best)
            }
            ModificationTag::FormulaTag(formula) => self.formula(formula, residue, term),
            ModificationTag::GlycanComposition(_)
            | ModificationTag::InfoTag(_)
            | ModificationTag::PositionConstraint(_) => Ok(None),
        }
    }
    pub(super) fn existing_full_id(
        &mut self,
        full_id: &str,
    ) -> Result<Option<Arc<ResidueModification>>> {
        self.lookup_charge(full_id)?;
        Ok(select_exact(self.database(), full_id, None, None, false).1)
    }
    pub(super) fn formula(
        &mut self,
        tag: &FormulaTag,
        residue: Option<char>,
        term: Option<TermSpecificity>,
    ) -> Result<Option<Arc<ResidueModification>>> {
        if tag.charge.is_some_and(|charge| charge != 0) {
            return Ok(None);
        }
        self.budget.text(&tag.formula_string)?;
        self.budget.consume(
            tag.formula_string
                .len()
                .saturating_add(1)
                .saturating_mul(256),
        )?;
        self.budget.allocate(
            tag.formula_string
                .len()
                .saturating_add(1)
                .saturating_mul(512)
                .saturating_add(2048),
        )?;
        // Parsing errors mean unresolved in source; resource checks above are
        // separate and must never be swallowed by this fallback.
        let formula = match EmpiricalFormula::parse(&tag.formula_string) {
            Ok(formula) => formula,
            Err(_) => return Ok(None),
        };
        if formula.is_empty() || formula.charge() != 0 {
            return Ok(None);
        }
        let specificity = term.unwrap_or(TermSpecificity::Anywhere);
        let (site, base) = match specificity {
            TermSpecificity::NTerm => (".n".to_owned(), None),
            TermSpecificity::CTerm => (".c".to_owned(), None),
            // These aren't reachable through the public resolution traversal.
            TermSpecificity::ProteinNTerm | TermSpecificity::ProteinCTerm => return Ok(None),
            TermSpecificity::Anywhere => {
                let Some(residue) = residue else {
                    return Ok(None);
                };
                if !residue.is_ascii() {
                    return Ok(None);
                }
                let base = if let Some(composition) = residue_composition(residue as u8) {
                    composition_formula(composition)
                        .checked_add(&composition_formula([0, 2, 0, 1, 0, 0]))?
                } else if matches!(residue, 'B' | 'Z' | 'X') {
                    EmpiricalFormula::default()
                } else {
                    return Ok(None);
                };
                (residue.to_string(), Some(base))
            }
        };
        self.budget.consume(
            formula
                .stored_atom_types()
                .saturating_add(12)
                .saturating_mul(256),
        )?;
        let canonical = canonical_formula(&formula, self.budget)?;
        self.budget.allocate(canonical.len().saturating_add(32))?;
        let full_id = format!("{site}[Formula:{canonical}]");
        self.lookup_charge(&full_id)?;
        if let (_, Some(existing), _) = select_exact(self.database(), &full_id, None, None, false) {
            return Ok(Some(existing));
        }
        let mono = formula.mono_mass();
        let average = formula.average_mass();
        let (absolute_mono, absolute_average) = match base {
            Some(base) => (mono + base.mono_mass(), average + base.average_mass()),
            None => (
                mono + composition_formula(if specificity == TermSpecificity::NTerm {
                    [0, 1, 0, 0, 0, 0]
                } else {
                    [0, 1, 0, 1, 0, 0]
                })
                .mono_mass(),
                0.,
            ),
        };
        self.budget.consume(full_id.len().saturating_add(2048))?;
        self.budget.allocate(full_id.len().saturating_add(2048))?;
        let magnitude = crate::param::value::format_float(mono.abs(), true);
        let record = ResidueModification::from_record(ModificationRecord {
            provenance: ModificationProvenance::MassOnly,
            full_id: full_id.clone(),
            full_name: format!("[{}{magnitude}]", if mono < 0. { "-" } else { "+" }),
            origin: if specificity == TermSpecificity::Anywhere {
                residue
            } else {
                None
            },
            term_specificity: specificity,
            diff_formula: formula,
            diff_mono_mass: mono,
            diff_average_mass: average,
            mono_mass: absolute_mono,
            average_mass: absolute_average,
            ..ModificationRecord::default()
        })?;
        self.insert_record(record)
    }
    /// Caller supplies an absent full ID; the shared transaction owns insertion.
    pub(super) fn insert_record(
        &mut self,
        record: ResidueModification,
    ) -> Result<Option<Arc<ResidueModification>>> {
        self.budget.record(&record)?;
        self.budget.allocate(
            record
                .payload_bytes()?
                .saturating_add(record.full_id().len())
                .saturating_add(1024),
        )?;
        if self.staged.is_none() {
            self.staged = Some(
                self.original
                    .clone_with_budget(&mut self.budget.work, &mut self.budget.bytes)?,
            );
        }
        let staged = self.staged.as_mut().expect("staged registry");
        self.budget.consume(staged.len())?;
        for existing in staged.entries() {
            self.budget.record(existing)?;
        }
        staged.charge_extension(&mut self.budget.work, &mut self.budget.bytes)?;
        staged.extend_records(vec![record])?;
        // FullId was absent and the append is atomic. Select the exact new Arc
        // without invoking the generic registry's ambiguity correction.
        Ok(staged.entries().last().cloned())
    }
}

fn source_matches(
    record: &ResidueModification,
    residue: Option<char>,
    term: Option<TermSpecificity>,
) -> bool {
    if term.is_some_and(|term| term != record.term_specificity()) {
        return false;
    }
    let query = residue.unwrap_or('?');
    let origin = record.origin().unwrap_or('X');
    if origin == 'X' {
        !(record.name().is_empty() && query != '?' && query != 'X')
    } else {
        origin == query || matches!(query, 'X' | '.' | '?')
    }
}

fn select_exact(
    registry: &ModificationsDB,
    name: &str,
    residue: Option<char>,
    term: Option<TermSpecificity>,
    defined: bool,
) -> (bool, Option<Arc<ResidueModification>>, usize) {
    let mut key_exists = false;
    let mut first = None;
    let mut count = 0;
    let mut consider = |record: &Arc<ResidueModification>| {
        key_exists = true;
        if source_matches(record, residue, term)
            && (!defined || record.provenance() == ModificationProvenance::Defined)
        {
            count += 1;
            if first.is_none() {
                first = Some(Arc::clone(record));
            }
        }
    };
    if name.is_empty() {
        // Source indexes empty short/full names and empty UniMod accessions.
        // The ordinary native registry intentionally omits those aliases.
        for record in registry.entries() {
            if record.name().is_empty()
                || record.full_name().is_empty()
                || record.record_id().is_none_or(|id| id == 0)
            {
                consider(record);
            }
        }
    } else if let Some(indices) = registry.exact_name_indices(name) {
        for &index in indices {
            consider(&registry.entries()[index]);
        }
    }
    (key_exists, first, count)
}

pub(super) fn canonical_formula(formula: &EmpiricalFormula, budget: &mut Budget) -> Result<String> {
    let count = formula.stored_atom_types();
    budget.consume(count.saturating_mul(1024))?;
    budget.allocate(count.saturating_mul(128).saturating_add(64))?;
    let mut elements = Vec::new();
    elements.try_reserve_exact(count).map_err(|_| limit())?;
    for (atom, &count) in &formula.atoms {
        let symbol = match atom.isotope {
            Some(isotope) => format!("({isotope}){}", atom.symbol),
            None => atom.symbol.to_owned(),
        };
        elements.push((symbol, count));
    }
    elements.sort_unstable_by(|a, b| a.0.cmp(&b.0));
    let mut text = String::new();
    text.try_reserve_exact(count.saturating_mul(32))
        .map_err(|_| limit())?;
    for (symbol, count) in elements {
        use std::fmt::Write;
        write!(text, "{symbol}{count}").map_err(|_| limit())?;
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_late_cumulative_limit_preserves_both_caller_objects() {
        let original = Peptidoform::parse("M[Formula:O]A[Formula:H2][not-found]").unwrap();
        let mut complete = original.clone();
        let mut complete_db = ModificationsDB::default();
        let mut budget = Budget::default();
        let warnings = resolve(&mut complete, &mut complete_db, &mut budget).unwrap();
        assert_eq!(complete_db.len(), 2);
        assert_eq!(warnings.len(), 1);
        let used = [
            MAX_PROFORMA_RESOLUTION_WORK - budget.work,
            MAX_PROFORMA_RESOLUTION_BYTES - budget.bytes,
            MAX_PROFORMA_RESOLUTION_ITEMS - budget.items,
        ];
        for (index, cost) in used.into_iter().enumerate() {
            let mut budget = Budget::default();
            match index {
                0 => budget.work = cost - 1,
                1 => budget.bytes = cost - 1,
                _ => budget.items = cost - 1,
            }
            let mut pf = original.clone();
            let mut db = ModificationsDB::default();
            assert!(resolve(&mut pf, &mut db, &mut budget).is_err());
            assert_eq!(pf, original);
            assert!(db.is_empty());
        }
    }

    #[test]
    fn formula_parse_failure_does_not_swallow_a_prior_resource_failure() {
        let mut pf = Peptidoform::parse("M[Formula:nonsense]").unwrap();
        let old = pf.clone();
        let mut db = ModificationsDB::default();
        let mut budget = Budget {
            work: 16,
            ..Budget::default()
        };
        assert!(resolve(&mut pf, &mut db, &mut budget).is_err());
        assert_eq!(pf, old);
        resolve(&mut pf, &mut db, &mut Budget::default()).unwrap();
        let SequenceSection::Element(element) = &pf.sequence[0] else {
            panic!()
        };
        assert!(element.modifications[0].resolved_mod.is_none());
        assert!(db.is_empty());
    }

    #[test]
    fn absent_queries_consume_work_even_in_an_empty_registry() {
        let db = ModificationsDB::default();
        let mut budget = Budget::default();
        let mut resolver = Resolver {
            original: &db,
            staged: None,
            budget: &mut budget,
            warnings: Vec::new(),
        };
        assert!(resolver.lookup("missing", None, None).unwrap().0.is_none());
        let spent = MAX_PROFORMA_RESOLUTION_WORK - resolver.budget.work;
        assert!(spent > 0);
        resolver.budget.work = spent - 1;
        assert!(resolver.lookup("missing", None, None).is_err());
        assert!(resolver.staged.is_none());
    }

    #[test]
    fn empty_operation_needs_no_copy_and_ignored_range_fields_are_not_read() {
        let mut db = ModificationsDB::global().clone();
        let mut empty = Peptidoform::default();
        resolve(
            &mut empty,
            &mut db,
            &mut Budget {
                work: 0,
                bytes: 0,
                items: 0,
            },
        )
        .unwrap();
        let element = SequenceElement {
            amino_acid: 'M',
            modifications: vec![Modification {
                alternatives: vec![(
                    ModificationTag::MassDelta(MassDelta {
                        mass: f64::NAN,
                        ..MassDelta::default()
                    }),
                    None,
                )],
                resolved_mod: None,
            }],
        };
        let mut pf = Peptidoform {
            sequence: vec![SequenceSection::ModifiedRange(ModifiedRange {
                elements: vec![element; 100],
                modifications: vec![],
            })],
            ..Peptidoform::default()
        };
        resolve(
            &mut pf,
            &mut db,
            &mut Budget {
                work: 1,
                bytes: 0,
                items: 1,
            },
        )
        .unwrap();
    }
}
