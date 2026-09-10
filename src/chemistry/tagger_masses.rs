// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Private Tagger mass dictionary from OpenMS4-core 7c029e8 Tagger.cpp and
//! Residue.cpp. Unlike peptide chemistry, this scalar dictionary retains the
//! source's finite signed masses, including B/Z/X's zero free-residue mass.
//! Ambiguous lookup uses first provider order, replacing allocation-dependent
//! source pointer order, without changing the registry's public lookup policy.
//! Empty short IDs resolve only actual anonymous records in a bounded private
//! scan; unrelated empty full-name/accession fields do not create aliases here.

use crate::chemistry::{
    EmpiricalFormula, ModificationsDB, ResidueModification, TermSpecificity, composition_formula,
    residue_composition,
};
use crate::{Error, Result};
use std::cmp::Ordering;

const NATURAL_19: &[u8; 19] = b"ACDEFGHKLMNPQRSTVWY";
const MAX_MODIFICATIONS: usize = 10_000;
const MAX_NAME_BYTES: usize = 1 << 20;
const MAX_REGISTRY_RECORDS: usize = 200_000;
const MAX_WORK: usize = 50_000_000;

/// Resolve owned numeric values only; no registry reference outlives construction.
/// Limits apply before allocating the result or collecting registry matches.
/// At most 10,019 numeric entries and 200,000 temporary borrowed match pointers
/// are retained. Both lookup passes share 1 MiB of query bytes and 50M charged
/// lookup, key-comparison, formula and table operations.
pub(super) fn resolve_mass_table(
    fixed: &[String],
    variable: &[String],
    registry: &ModificationsDB,
) -> Result<Vec<(f64, u8)>> {
    let mut work = MAX_WORK;
    resolve_with_work(fixed, variable, registry, &mut work)
}

fn resolve_with_work(
    fixed: &[String],
    variable: &[String],
    registry: &ModificationsDB,
    work: &mut usize,
) -> Result<Vec<(f64, u8)>> {
    let count = fixed
        .len()
        .checked_add(variable.len())
        .filter(|&n| n <= MAX_MODIFICATIONS)
        .ok_or_else(|| invalid("modification count limit exceeded"))?;
    if count != 0 && registry.len() > MAX_REGISTRY_RECORDS {
        return Err(invalid("registry record limit exceeded"));
    }
    let mut name_bytes = 0usize;
    for name in fixed.iter().chain(variable) {
        add_name_bytes(&mut name_bytes, name)?;
    }
    charge(work, 19 * 4096 + count)?;
    let water = composition_formula([0, 2, 0, 1, 0, 0]);
    let water_mass = water.mono_mass();
    let mut table = Vec::new();
    table
        .try_reserve_exact(NATURAL_19.len() + count)
        .map_err(|_| invalid("mass table allocation failed"))?;
    for &letter in NATURAL_19 {
        let mass = free_formula(letter, &water)?.mono_mass() - water_mass;
        insert(&mut table, mass, letter);
    }
    for (index, name) in fixed.iter().chain(variable).enumerate() {
        let initial = lookup(registry, name, None, None, work)?;
        let origin = initial.origin().unwrap_or('X');
        let formula = free_formula(origin as u8, &water)?;
        // Residue::setModification(short_id) performs a second lookup. Even a
        // full-ID request must not bypass it or carry its first chemistry through.
        add_name_bytes(&mut name_bytes, initial.name())?;
        let applied = lookup(
            registry,
            initial.name(),
            Some(origin),
            Some(TermSpecificity::Anywhere),
            work,
        )?;
        let mass = modified_internal_mass(formula, applied, water_mass, work)?;
        // Covers first-value search, removal shift, insertion shift and binary
        // comparisons, including repeated fixed assignments to one letter.
        charge(work, table.len().saturating_mul(4).saturating_add(1))?;
        if index < fixed.len() {
            if let Some(position) = table.iter().position(|&(_, letter)| letter == origin as u8) {
                table.remove(position);
            }
        }
        insert(&mut table, mass, origin as u8);
    }
    Ok(table)
}

fn add_name_bytes(total: &mut usize, name: &str) -> Result<()> {
    *total = total
        .checked_add(name.len())
        .filter(|&n| n <= MAX_NAME_BYTES)
        .ok_or_else(|| invalid("modification query byte limit exceeded"))?;
    Ok(())
}

fn lookup<'a>(
    registry: &'a ModificationsDB,
    name: &str,
    origin: Option<char>,
    term: Option<TermSpecificity>,
    work: &mut usize,
) -> Result<&'a ResidueModification> {
    // A BTree lookup needs fewer than 1024 string comparisons on a usize-sized
    // map. Charge query bytes per comparison, plus all possible match visits,
    // before the registry allocates its borrowed-match vector. Long unrelated
    // keys cost at most the query length in each comparison.
    charge(
        work,
        name.len()
            .saturating_add(1)
            .saturating_mul(1024)
            .saturating_add(registry.len().saturating_mul(2)),
    )?;
    let matches_origin = |record: &&ResidueModification| {
        origin.is_none_or(|query| {
            let target = record.origin().unwrap_or('X');
            if target == 'X' {
                // Source anonymous X means an actual X, not any residue.
                !record.name().is_empty() || query == 'X'
            } else {
                target == query || query == 'X'
            }
        })
    };
    let found = if name.is_empty() {
        // The ordinary native registry intentionally omits empty aliases. Source
        // anonymous records nevertheless require a second empty-short-ID lookup.
        registry
            .entries()
            .iter()
            .map(|record| record.as_ref())
            .filter(|record| {
                record.name().is_empty()
                    && term.is_none_or(|term| record.term_specificity() == term)
            })
            .find(matches_origin)
    } else {
        registry
            .find(name, None, term)
            .into_iter()
            .find(matches_origin)
    };
    found.ok_or_else(|| invalid("modification does not resolve for the requested residue/terminus"))
}

fn free_formula(letter: u8, water: &EmpiricalFormula) -> Result<EmpiricalFormula> {
    match residue_composition(letter) {
        Some(composition) => composition_formula(composition).checked_add(water),
        None if matches!(letter, b'B' | b'Z' | b'X') => Ok(EmpiricalFormula::default()),
        // Source's byte-indexed ResidueDB returns null for an unknown origin.
        None => Err(invalid("modification origin is not a source residue")),
    }
}

fn modified_internal_mass(
    formula: EmpiricalFormula,
    modification: &ResidueModification,
    water_mass: f64,
    work: &mut usize,
) -> Result<f64> {
    let delta = modification.diff_formula();
    let absolute = modification.absolute_formula();
    let cells = formula
        .atoms
        .len()
        .saturating_add(delta.atoms.len())
        .saturating_add(absolute.map_or(0, |f| f.atoms.len()));
    // Includes formula-map copying/combining, element lookup and isotope scans.
    charge(work, cells.saturating_add(1).saturating_mul(256))?;
    let mut mass = formula.mono_mass();
    let changes = modification.diff_mono_mass() != 0.0
        || modification.diff_average_mass() != 0.0
        || !delta.is_empty();
    if changes {
        if !delta.is_empty() {
            mass = formula.checked_add(delta)?.mono_mass();
        } else if let Some(absolute) = absolute {
            mass = absolute.mono_mass();
        } else if modification.mono_mass() != 0.0 {
            mass = modification.mono_mass();
        } else if modification.diff_mono_mass() != 0.0 {
            mass += modification.diff_mono_mass();
        }
    }
    // Exact source order is free-residue mass followed by water subtraction.
    // Finite zero/negative values are valid private dictionary keys.
    let mass = mass - water_mass;
    if mass.is_finite() {
        Ok(mass)
    } else {
        Err(invalid("calculated residue mass is not finite"))
    }
}

fn insert(table: &mut Vec<(f64, u8)>, mass: f64, letter: u8) {
    match table.binary_search_by(|&(stored, _)| {
        if stored == mass {
            Ordering::Equal // std::map treats positive/negative zero as one key.
        } else {
            stored.total_cmp(&mass)
        }
    }) {
        Ok(index) => table[index].1 = letter,
        Err(index) => table.insert(index, (mass, letter)),
    }
}

fn charge(work: &mut usize, amount: usize) -> Result<()> {
    *work = work
        .checked_sub(amount)
        .ok_or_else(|| invalid("construction work limit exceeded"))?;
    Ok(())
}

fn invalid(message: &str) -> Error {
    Error::InvalidValue(format!("Tagger {message}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::{AASequence, ModificationRecord};

    const H: f64 = 1.007_825_031_9;
    const N: f64 = 14.003_074;
    const O: f64 = 15.994_915;
    const S: f64 = 31.972_070_73;

    fn free_mass(counts: [u32; 5]) -> f64 {
        [12.0, H, N, O, S]
            .into_iter()
            .zip(counts)
            .fold(0.0, |sum, (mass, count)| sum + mass * f64::from(count))
    }
    fn water_mass() -> f64 {
        H * 2.0 + O
    }
    fn record(name: &str, origin: char, delta: f64) -> ModificationRecord {
        ModificationRecord {
            name: name.into(),
            origin: Some(origin),
            diff_mono_mass: delta,
            ..Default::default()
        }
    }
    fn database(records: Vec<ModificationRecord>) -> ModificationsDB {
        ModificationsDB::from_records(
            records
                .into_iter()
                .map(|r| ResidueModification::from_record(r).unwrap())
                .collect(),
        )
        .unwrap()
    }
    fn table(fixed: &[&str], variable: &[&str], db: &ModificationsDB) -> Vec<(f64, u8)> {
        resolve_mass_table(
            &fixed.iter().map(|s| (*s).into()).collect::<Vec<_>>(),
            &variable.iter().map(|s| (*s).into()).collect::<Vec<_>>(),
            db,
        )
        .unwrap()
    }
    fn masses(table: &[(f64, u8)], letter: u8) -> Vec<f64> {
        table
            .iter()
            .filter_map(|&(mass, code)| (code == letter).then_some(mass))
            .collect()
    }

    #[test]
    fn nineteen_source_residues_use_free_formula_then_remove_water() {
        // Full residue counts transcribed independently from ResidueDB.cpp.
        let counts = [
            [3, 7, 1, 2, 0],
            [3, 7, 1, 2, 1],
            [4, 7, 1, 4, 0],
            [5, 9, 1, 4, 0],
            [9, 11, 1, 2, 0],
            [2, 5, 1, 2, 0],
            [6, 9, 3, 2, 0],
            [6, 14, 2, 2, 0],
            [6, 13, 1, 2, 0],
            [5, 11, 1, 2, 1],
            [4, 8, 2, 3, 0],
            [5, 9, 1, 2, 0],
            [5, 10, 2, 3, 0],
            [6, 14, 4, 2, 0],
            [3, 7, 1, 3, 0],
            [4, 9, 1, 3, 0],
            [5, 11, 1, 2, 0],
            [11, 12, 2, 2, 0],
            [9, 11, 1, 3, 0],
        ];
        let result = table(&[], &[], &ModificationsDB::default());
        assert_eq!(result.len(), 19);
        assert!(result.windows(2).all(|pair| pair[0].0 < pair[1].0));
        for (&letter, counts) in NATURAL_19.iter().zip(counts) {
            assert_eq!(masses(&result, letter), [free_mass(counts) - water_mass()]);
        }
        assert!(masses(&result, b'I').is_empty());
        let alanine = AASequence::parse("A")
            .unwrap()
            .internal_mass(0)
            .unwrap()
            .unwrap();
        assert_ne!(masses(&result, b'A')[0].to_bits(), alanine.to_bits());
    }

    #[test]
    fn named_modifications_use_formula_mass_and_fixed_then_variable_order() {
        let result = table(
            &["Carbamidomethyl (C)"],
            &["Oxidation (M)"],
            ModificationsDB::global(),
        );
        assert_eq!(
            masses(&result, b'C'),
            [free_mass([5, 10, 2, 3, 1]) - water_mass()]
        );
        assert_eq!(
            masses(&result, b'M'),
            [
                free_mass([5, 11, 1, 2, 1]) - water_mass(),
                free_mass([5, 11, 1, 3, 1]) - water_mass(),
            ]
        );
    }

    #[test]
    fn fixed_assignments_restart_from_base_and_variables_follow_all_fixed_changes() {
        let db = database(vec![record("Low", 'M', -5.0), record("High", 'M', 12.5)]);
        let result = table(&["Low", "High"], &["Low"], &db);
        let base = free_mass([5, 11, 1, 2, 1]);
        assert_eq!(
            masses(&result, b'M'),
            [base - 5.0 - water_mass(), base + 12.5 - water_mass()]
        );
    }

    #[test]
    fn fixed_i_does_not_remove_l_and_exact_collisions_overwrite_last() {
        let db = database(vec![record("ExtraI", 'I', 1.0), record("PlainI", 'I', 0.0)]);
        let result = table(&["ExtraI"], &[], &db);
        assert_eq!(masses(&result, b'L').len(), 1);
        assert_eq!(masses(&result, b'I').len(), 1);
        let result = table(&["PlainI"], &[], &db);
        assert!(masses(&result, b'L').is_empty());
        assert_eq!(masses(&result, b'I').len(), 1);
        let mut result = Vec::new();
        insert(&mut result, -0.0, b'A');
        insert(&mut result, 0.0, b'C');
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].1, b'C');
        assert_eq!(result[0].0.to_bits(), (-0.0f64).to_bits());
    }

    #[test]
    fn same_short_id_is_resolved_again_and_provider_order_is_deterministic() {
        let mut first = record("Lab", 'M', 1.0);
        first.full_id = "Lab first".into();
        let mut second = record("Lab", 'M', 9.0);
        second.full_id = "Lab second".into();
        let db = database(vec![first.clone(), second.clone()]);
        assert!(db.get_modification("Lab", None, None).is_err());
        let result = table(&["Lab second"], &[], &db);
        let base = free_mass([5, 11, 1, 2, 1]);
        assert_eq!(masses(&result, b'M'), [base + 1.0 - water_mass()]);
        let reversed = database(vec![second.clone(), first.clone()]);
        assert_eq!(
            masses(&table(&["Lab first"], &[], &reversed), b'M'),
            [base + 9.0 - water_mass()]
        );
        // Identical full IDs do not make distinct caller chemistry identical.
        // The first lookup reaches the second record by its distinct full name;
        // the required short-ID lookup then deliberately chooses the first.
        first.full_id = "Same ID".into();
        second.full_id = "Same ID".into();
        second.full_name = "Second chemical record".into();
        let db = database(vec![first, second]);
        assert_eq!(
            masses(&table(&["Second chemical record"], &[], &db), b'M'),
            [base + 1.0 - water_mass()]
        );
    }

    #[test]
    fn terminal_request_can_resolve_to_anywhere_but_keeps_original_letter() {
        let mut terminal = record("Shared", 'M', 99.0);
        terminal.term_specificity = TermSpecificity::NTerm;
        let db = database(vec![terminal.clone(), record("Shared", 'M', 2.0)]);
        let base = free_mass([5, 11, 1, 2, 1]);
        assert_eq!(
            masses(&table(&["Shared (N-term M)"], &[], &db), b'M'),
            [base + 2.0 - water_mass()]
        );
        let only_terminal = database(vec![terminal]);
        assert!(resolve_mass_table(&["Shared (N-term M)".into()], &[], &only_terminal).is_err());
        let mut wildcard = record("Shared", 'X', 99.0);
        wildcard.term_specificity = TermSpecificity::NTerm;
        let db = database(vec![wildcard, record("Shared", 'M', 25.0)]);
        assert_eq!(
            masses(&table(&["Shared (N-term)"], &[], &db), b'X'),
            [25.0 - water_mass()]
        );
    }

    #[test]
    fn custom_mass_and_formula_precedence_matches_free_residue_source_arithmetic() {
        let base = free_mass([5, 11, 1, 2, 1]);
        let mut absolute = record("Absolute", 'M', 1.0);
        absolute.mono_mass = 300.0;
        let mut no_op = record("NoOp", 'M', 0.0);
        no_op.mono_mass = 300.0;
        no_op.absolute_formula = Some("C2H5NO2".parse().unwrap());
        let mut formula = absolute.clone();
        formula.name = "Formula".into();
        formula.diff_formula = "O".parse().unwrap();
        let mut full = absolute.clone();
        full.name = "Full".into();
        full.absolute_formula = Some("C2H5NO2".parse().unwrap());
        let mut average_only = record("Average", 'M', 0.0);
        average_only.diff_average_mass = 2.0;
        let db = database(vec![
            record("Delta", 'M', 12.5),
            absolute,
            no_op,
            formula,
            full,
            average_only,
        ]);
        for (name, mass) in [
            ("Delta", base + 12.5),
            ("Absolute", 300.0),
            ("NoOp", base),
            ("Formula", free_mass([5, 11, 1, 3, 1])),
            ("Full", free_mass([2, 5, 1, 2, 0])),
            ("Average", base),
        ] {
            assert_eq!(
                masses(&table(&[name], &[], &db), b'M'),
                [mass - water_mass()],
                "{name}"
            );
        }
    }

    #[test]
    fn ambiguous_letters_and_signed_custom_values_keep_source_scalar_semantics() {
        let mut records = Vec::new();
        for letter in ['B', 'Z', 'X'] {
            records.push(record(&letter.to_string(), letter, 0.0));
        }
        let mut negative = record("Negative", 'A', 1.0);
        negative.mono_mass = -100.0;
        let mut zero = record("Zero", 'C', 1.0);
        zero.mono_mass = water_mass();
        records.extend([negative, zero]);
        let db = database(records);
        for letter in ['B', 'Z', 'X'] {
            assert_eq!(
                masses(&table(&[], &[&letter.to_string()], &db), letter as u8),
                [-water_mass()]
            );
        }
        let result = table(&["Negative", "Zero"], &[], &db);
        assert_eq!(masses(&result, b'A'), [-100.0 - water_mass()]);
        assert_eq!(masses(&result, b'C'), [0.0]);
    }

    #[test]
    fn charge_only_noop_formula_and_explicit_empty_absolute_are_distinct() {
        let mut noop = record("Charge", 'A', 0.0);
        noop.diff_formula = "+".parse().unwrap();
        let mut empty = record("Empty", 'A', 1.0);
        empty.absolute_formula = Some(EmpiricalFormula::default());
        let db = database(vec![noop, empty]);
        assert_eq!(
            masses(&table(&["Charge"], &[], &db), b'A'),
            [free_mass([3, 7, 1, 2, 0]) - water_mass()]
        );
        assert_eq!(masses(&table(&["Empty"], &[], &db), b'A'), [-water_mass()]);
    }

    #[test]
    fn unknown_origins_missing_names_and_budget_fail_before_partial_output() {
        let db = database(vec![record("Valid", 'V', 0.0)]);
        assert!(resolve_mass_table(&["Missing".into()], &[], &db).is_err());
        let unchanged = db.entries().to_vec();
        assert!(resolve_with_work(&["Valid".into()], &[], &db, &mut 0).is_err());
        assert_eq!(db.entries(), unchanged);
        assert!(ResidueModification::from_record(record("Bad origin", '?', 0.0)).is_err());
        assert!(free_formula(b'?', &"H2O".parse().unwrap()).is_err());
        let long = "x".repeat(MAX_NAME_BYTES + 1);
        assert!(resolve_mass_table(&[long], &[], &db).is_err());
        assert!(resolve_mass_table(&vec![String::new(); MAX_MODIFICATIONS + 1], &[], &db).is_err());
        let mut work = 50;
        assert!(lookup(&db, "Valid", None, None, &mut work).is_err());
        assert_eq!(work, 50);
        let mut overflow = record("Overflow", 'A', 1.0);
        overflow.diff_formula = "C2147483647".parse().unwrap();
        let db = database(vec![overflow]);
        assert!(resolve_mass_table(&["Overflow".into()], &[], &db).is_err());
    }

    #[test]
    fn variable_exact_mass_collision_keeps_last_requested_letter() {
        let mut alanine = record("LabA", 'A', 1.0);
        alanine.mono_mass = 300.0;
        let mut cysteine = record("LabC", 'C', 1.0);
        cysteine.mono_mass = 300.0;
        let db = database(vec![alanine, cysteine]);
        let result = table(&[], &["LabA", "LabC"], &db);
        assert_eq!(result.last(), Some(&(300.0 - water_mass(), b'C')));
        let result = table(&[], &["LabC", "LabA"], &db);
        assert_eq!(result.last(), Some(&(300.0 - water_mass(), b'A')));
    }

    #[test]
    fn anonymous_short_ids_reresolve_privately_without_global_empty_aliases() {
        let mut unknown = record("", 'X', 1.0);
        unknown.full_id = "X[100]".into();
        unknown.mono_mass = 100.0 + water_mass();
        let mut first = record("", 'M', 1.0);
        first.full_id = "M[+1]".into();
        let mut second = record("", 'M', 5.0);
        second.full_id = "M[+5]".into();
        let db = database(vec![record("Named", 'M', 99.0), unknown, first, second]);
        assert!(db.find("", None, None).is_empty());
        let result = table(&["M[+5]"], &["X[100]"], &db);
        // The anonymous X entry does not wildcard-match an M query; the named
        // record's empty full name does not qualify as an empty short ID.
        assert_eq!(
            masses(&result, b'M'),
            [free_mass([5, 11, 1, 2, 1]) + 1.0 - water_mass()]
        );
        assert_eq!(masses(&result, b'X'), [100.0]);
        drop(db);
        assert_eq!(masses(&result, b'X'), [100.0]);
    }
}
