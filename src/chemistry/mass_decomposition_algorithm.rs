// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded source MassDecompositionAlgorithm with private integer residue tables.
//! The alphabet is ordered by symbol, and the real-mass search retains the
//! source's exclusive upper integer bound. See MASS_DECOMPOSITION_ALGORITHM_SUPPORT.md.

use super::{
    MassDecomposition, ModificationsDB, ResidueModification, composition_formula,
    residue_composition,
};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::str::FromStr;

pub const MAX_DECOMPOSITION_WORK: usize = 50_000_000;
pub const MAX_DECOMPOSITION_TABLE_CELLS: usize = 4_000_000;
pub const MAX_DECOMPOSITION_OUTPUTS: usize = 100_000;
pub const MAX_DECOMPOSITION_BYTES: usize = 64 * 1024 * 1024;
const MAX_NAMES: usize = 10_000;
const MAX_NAME_BYTES: usize = 1024 * 1024;
const MAX_REGISTRY_RECORDS: usize = 200_000;
const MAX_INTEGER: u64 = (1 << 53) - 1;

/// All eight named sets from source ResidueDB. Natural19J actually contains the
/// same twenty residues as Natural20; the source does not insert J into it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DecompositionResidueSet {
    All,
    AllNatural,
    Ambiguous,
    AmbiguousWithoutX,
    Natural19J,
    #[default]
    Natural19WithoutI,
    Natural19WithoutL,
    Natural20,
}
impl DecompositionResidueSet {
    pub const ALL: [Self; 8] = [
        Self::All,
        Self::AllNatural,
        Self::Ambiguous,
        Self::AmbiguousWithoutX,
        Self::Natural19J,
        Self::Natural19WithoutI,
        Self::Natural19WithoutL,
        Self::Natural20,
    ];
    pub fn name(self) -> &'static str {
        match self {
            Self::All => "All",
            Self::AllNatural => "AllNatural",
            Self::Ambiguous => "Ambiguous",
            Self::AmbiguousWithoutX => "AmbiguousWithoutX",
            Self::Natural19J => "Natural19J",
            Self::Natural19WithoutI => "Natural19WithoutI",
            Self::Natural19WithoutL => "Natural19WithoutL",
            Self::Natural20 => "Natural20",
        }
    }
    fn symbols(self) -> &'static [u8] {
        match self {
            Self::All | Self::Ambiguous => b"ABCDEFGHIJKLMNOPQRSTUVWXYZ",
            Self::AmbiguousWithoutX => b"ABCDEFGHIJKLMNOPQRSTUVWYZ",
            Self::AllNatural => b"ACDEFGHIKLMNOPQRSTUVWY",
            Self::Natural19J | Self::Natural20 => b"ACDEFGHIKLMNPQRSTVWY",
            Self::Natural19WithoutI => b"ACDEFGHKLMNPQRSTVWY",
            Self::Natural19WithoutL => b"ACDEFGHIKMNPQRSTVWY",
        }
    }
}
impl FromStr for DecompositionResidueSet {
    type Err = Error;
    fn from_str(value: &str) -> Result<Self> {
        Self::ALL
            .into_iter()
            .find(|v| v.name() == value)
            .ok_or_else(|| invalid("unknown residue set"))
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MassDecompositionOptions {
    pub decomp_weights_precision: f64,
    pub tolerance: f64,
    pub fixed_modifications: Vec<String>,
    pub variable_modifications: Vec<String>,
    pub residue_set: DecompositionResidueSet,
}
impl Default for MassDecompositionOptions {
    fn default() -> Self {
        Self {
            decomp_weights_precision: 0.01,
            tolerance: 0.3,
            fixed_modifications: Vec::new(),
            variable_modifications: Vec::new(),
            residue_set: DecompositionResidueSet::default(),
        }
    }
}

/// Owns resolved numeric masses, so a custom registry may be dropped afterward.
/// Construction and each query have separate 50M work / 64MiB logical allocation
/// budgets. Allocator metadata, fragmentation and caller-reserved capacity are
/// outside the conservative payload model; this is not a physical RSS ceiling.
#[derive(Debug)]
pub struct MassDecompositionAlgorithm {
    options: MassDecompositionOptions,
    alphabet: Vec<(char, f64)>,
    diagnostics: Vec<String>,
    solver: IntegerSolver,
    precision: f64,
    minimum_error: f64,
    maximum_error: f64,
}
impl MassDecompositionAlgorithm {
    pub fn new() -> Result<Self> {
        Self::from_options(MassDecompositionOptions::default())
    }
    pub fn from_options(options: MassDecompositionOptions) -> Result<Self> {
        Self::with_registry(options, ModificationsDB::global())
    }
    pub fn with_registry(
        options: MassDecompositionOptions,
        registry: &ModificationsDB,
    ) -> Result<Self> {
        if !options.decomp_weights_precision.is_finite()
            || options.decomp_weights_precision <= 0.0
            || !options.tolerance.is_finite()
            || options.tolerance < 0.0
        {
            return Err(invalid(
                "precision must be positive finite and tolerance nonnegative finite",
            ));
        }
        let mut budget = Budget::default();
        let mut diagnostics = Vec::new();
        let alphabet = alphabet(&options, registry, &mut budget, &mut diagnostics)?;
        budget.allocate(alphabet.len() * 8)?;
        let mut weights = Vec::with_capacity(alphabet.len());
        for &(_, mass) in &alphabet {
            if !mass.is_finite() || mass <= 0.0 {
                return Err(invalid(
                    "every resolved alphabet mass must be positive finite",
                ));
            }
            let integer = integer((mass / options.decomp_weights_precision + 0.5).floor())?;
            if integer == 0 {
                return Err(invalid("precision rounds an alphabet mass to zero"));
            }
            weights.push(integer);
        }
        // Do not recompute weights after GCD scaling: rounding could differ.
        let divisor = weights.iter().copied().reduce(gcd).unwrap_or(1);
        let precision = options.decomp_weights_precision * divisor as f64;
        if !precision.is_finite() || precision <= 0.0 {
            return Err(invalid("scaled precision is not finite positive"));
        }
        for weight in &mut weights {
            *weight /= divisor;
        }
        let mut minimum_error = 0.0f64;
        let mut maximum_error = 0.0f64;
        for (&(_, mass), &weight) in alphabet.iter().zip(&weights) {
            let error = (precision * weight as f64 - mass) / mass;
            if !error.is_finite() {
                return Err(invalid("rounding error is not finite"));
            }
            minimum_error = minimum_error.min(error);
            maximum_error = maximum_error.max(error);
        }
        let solver = IntegerSolver::new(weights, &mut budget)?;
        Ok(Self {
            options,
            alphabet,
            diagnostics,
            solver,
            precision,
            minimum_error,
            maximum_error,
        })
    }
    pub fn options(&self) -> &MassDecompositionOptions {
        &self.options
    }
    /// Alphabetical symbols, including the source lowercase variable labels.
    pub fn alphabet(&self) -> &[(char, f64)] {
        &self.alphabet
    }
    /// Ignored ambiguous/massless modifications, returned instead of stderr.
    pub fn diagnostics(&self) -> &[String] {
        &self.diagnostics
    }
    pub fn set_options(&mut self, options: MassDecompositionOptions) -> Result<()> {
        self.set_options_with_registry(options, ModificationsDB::global())
    }
    pub fn set_options_with_registry(
        &mut self,
        options: MassDecompositionOptions,
        registry: &ModificationsDB,
    ) -> Result<()> {
        let next = Self::with_registry(options, registry)?;
        *self = next;
        Ok(())
    }
    pub fn decompositions(&self, mass: f64) -> Result<Vec<MassDecomposition>> {
        let mut output = Vec::new();
        self.append_decompositions(&mut output, mass)?;
        Ok(output)
    }
    /// Append in source integer-mass/recursive order. Existing contents and their
    /// allocation are untouched on failure, including a late output/work limit.
    pub fn append_decompositions(
        &self,
        output: &mut Vec<MassDecomposition>,
        mass: f64,
    ) -> Result<()> {
        self.append_with_budget(output, mass, &mut Budget::default())
    }
    fn append_with_budget(
        &self,
        output: &mut Vec<MassDecomposition>,
        mass: f64,
        budget: &mut Budget,
    ) -> Result<()> {
        if !mass.is_finite() || mass < 0.0 {
            return Err(invalid("query mass must be finite nonnegative"));
        }
        if output.len() > MAX_DECOMPOSITION_OUTPUTS {
            return Err(invalid("output count limit exceeded"));
        }
        budget.charge(output.len())?;
        let mut output_bytes = 0usize;
        for value in output.iter() {
            output_bytes = output_bytes
                .checked_add(value.payload_bytes())
                .ok_or_else(|| invalid("output size overflow"))?;
        }
        if output_bytes > MAX_DECOMPOSITION_BYTES {
            return Err(invalid("existing output byte limit exceeded"));
        }
        let start = integer(
            ((1.0 + self.minimum_error) * (mass - self.options.tolerance) / self.precision).ceil(),
        )?;
        let end = integer(
            ((1.0 + self.maximum_error) * (mass + self.options.tolerance) / self.precision).floor(),
        )?;
        budget.charge(
            usize::try_from(end.saturating_sub(start))
                .map_err(|_| invalid("integer search interval too large"))?,
        )?;
        let mut found = Vec::new();
        budget.allocate(self.alphabet.len() * 4)?;
        let mut counts = vec![0u32; self.alphabet.len()];
        for integer_mass in start..end {
            self.solver.collect(
                integer_mass,
                self.alphabet.len() - 1,
                &mut counts,
                budget,
                &mut |counts, budget| {
                    budget.charge(counts.len())?;
                    let mut parent_mass = 0.0;
                    for (&(_, weight), &count) in self.alphabet.iter().zip(counts) {
                        parent_mass += weight * f64::from(count);
                    }
                    if !parent_mass.is_finite() {
                        return Err(invalid("reconstructed mass overflow"));
                    }
                    if (parent_mass - mass).abs() > self.options.tolerance {
                        return Ok(());
                    }
                    if output.len() + found.len() >= MAX_DECOMPOSITION_OUTPUTS {
                        return Err(invalid("output count limit exceeded"));
                    }
                    let keys = counts.iter().filter(|&&c| c != 0).count();
                    let bytes = std::mem::size_of::<MassDecomposition>()
                        + if keys == 0 { 0 } else { 256 + keys * 64 };
                    output_bytes = output_bytes
                        .checked_add(bytes)
                        .filter(|&b| b <= MAX_DECOMPOSITION_BYTES)
                        .ok_or_else(|| invalid("output byte limit exceeded"))?;
                    budget.allocate(bytes + 2 * std::mem::size_of::<MassDecomposition>())?;
                    budget.charge(counts.len() * 16 + 1)?;
                    found
                        .try_reserve(1)
                        .map_err(|_| invalid("output allocation failed"))?;
                    found.push(MassDecomposition::from_counts(
                        self.alphabet
                            .iter()
                            .zip(counts)
                            .map(|(&(symbol, _), &count)| (symbol as u8, count)),
                    )?);
                    Ok(())
                },
            )?;
        }
        // Reserve only after every fallible scientific/budget operation. A
        // failed try_reserve leaves the existing vector and values unchanged.
        if !found.is_empty() {
            budget.allocate(output.len() * std::mem::size_of::<MassDecomposition>())?;
        }
        output
            .try_reserve(found.len())
            .map_err(|_| invalid("append allocation failed"))?;
        output.append(&mut found);
        Ok(())
    }
}

fn alphabet(
    options: &MassDecompositionOptions,
    registry: &ModificationsDB,
    budget: &mut Budget,
    diagnostics: &mut Vec<String>,
) -> Result<Vec<(char, f64)>> {
    let name_count = options
        .fixed_modifications
        .len()
        .checked_add(options.variable_modifications.len())
        .filter(|&n| n <= MAX_NAMES)
        .ok_or_else(|| invalid("modification count limit exceeded"))?;
    budget.charge(name_count)?;
    let mut name_bytes = 0usize;
    for name in options
        .fixed_modifications
        .iter()
        .chain(&options.variable_modifications)
    {
        name_bytes = name_bytes
            .checked_add(name.len())
            .filter(|&n| n <= MAX_NAME_BYTES)
            .ok_or_else(|| invalid("modification name byte limit exceeded"))?;
    }
    if name_count != 0 && registry.len() > MAX_REGISTRY_RECORDS {
        return Err(invalid("registry record limit exceeded"));
    }
    budget.allocate(name_bytes + name_count * std::mem::size_of::<String>())?;
    budget.charge(26 * 4096)?;
    let water = composition_formula([0, 2, 0, 1, 0, 0]);
    budget.allocate(128 * 128)?; // bounded byte-symbol map, masses and diagnostics slots
    let mut table = BTreeMap::new();
    for &symbol in options.residue_set.symbols() {
        let mass = match residue_composition(symbol) {
            Some(formula) => composition_formula(formula)
                .checked_add(&water)?
                .mono_mass(),
            None => 0.0,
        } - water.mono_mass();
        table.insert(symbol, mass);
    }
    for (variable, names) in [
        (false, &options.fixed_modifications),
        (true, &options.variable_modifications),
    ] {
        budget.allocate(512)?;
        let mut definitions = BTreeMap::<&str, &ResidueModification>::new();
        for name in names {
            budget.charge(
                (name.len() + 1)
                    .saturating_mul(1024)
                    .saturating_add(registry.len().saturating_mul(2)),
            )?;
            budget.allocate(
                registry
                    .len()
                    .saturating_mul(2 * std::mem::size_of::<&ResidueModification>()),
            )?;
            let record = registry
                .find(name, None, None)
                .into_iter()
                .next()
                .ok_or_else(|| invalid("modification name does not resolve"))?;
            budget.charge((record.full_id().len() + 1).saturating_mul(definitions.len() + 1))?;
            budget.allocate(128)?;
            definitions.entry(record.full_id()).or_insert(record);
        }
        if variable && definitions.len() > 26 {
            return Err(invalid(
                "source supports at most26 distinct variable modifications",
            ));
        }
        for (index, record) in definitions.values().enumerate() {
            let symbol = if variable {
                b'a' + index as u8
            } else {
                record.origin().unwrap_or('X') as u8
            };
            let origin = record.origin().unwrap_or('X');
            if !origin.is_ascii() {
                return Err(invalid("modification origin must be an ASCII byte"));
            }
            let ignored = if origin == 'X' {
                Some("ambiguous origin")
            } else if record.mono_mass() == 0.0 && record.diff_mono_mass() == 0.0 {
                Some("no monoisotopic mass")
            } else {
                None
            };
            if let Some(reason) = ignored {
                budget.charge(record.full_id().len() + 128)?;
                budget.allocate(record.full_id().len() + 128)?;
                diagnostics.push(format!("Ignored {}: {reason}", record.full_id()));
                continue;
            }
            let mass = if record.mono_mass() != 0.0 {
                record.mono_mass()
            } else {
                *table.entry(origin as u8).or_insert(0.0) + record.diff_mono_mass()
            };
            table.insert(symbol, mass);
        }
    }
    budget.allocate(table.len() * 128)?;
    Ok(table
        .into_iter()
        .map(|(symbol, mass)| (char::from(symbol), mass))
        .collect())
}

#[derive(Debug)]
struct IntegerSolver {
    weights: Vec<u64>,
    table: Vec<Vec<u64>>,
    lcms: Vec<u64>,
    mass_in_lcms: Vec<u64>,
    infinity: u64,
}
impl IntegerSolver {
    fn new(weights: Vec<u64>, budget: &mut Budget) -> Result<Self> {
        if weights.is_empty() || weights.len() > 128 || weights.contains(&0) {
            return Err(invalid("invalid integer alphabet"));
        }
        let base = usize::try_from(weights[0]).map_err(|_| invalid("table width overflow"))?;
        let cells = base
            .checked_mul(weights.len())
            .filter(|&n| n <= MAX_DECOMPOSITION_TABLE_CELLS)
            .ok_or_else(|| invalid("precision-derived residue table limit exceeded"))?;
        budget.charge(cells)?;
        budget.allocate(cells * 8 + weights.len() * 64)?;
        let infinity = weights[0]
            .checked_mul(*weights.last().unwrap())
            .ok_or_else(|| invalid("source infinity sentinel overflow"))?;
        let mut table = vec![vec![infinity; base]; weights.len()];
        for row in &mut table {
            row[0] = 0;
        }
        let mut lcms = vec![0; weights.len()];
        let mut mass_in_lcms = vec![0; weights.len()];
        for i in 1..weights.len() {
            let divisor = gcd(weights[0], weights[i]);
            lcms[i] = weights[i]
                .checked_mul(weights[0])
                .ok_or_else(|| invalid("least-common-multiple product overflow"))?
                / divisor;
            mass_in_lcms[i] = weights[0] / divisor;
        }
        if weights.len() >= 2 {
            let increment = (weights[1] % weights[0]) as usize;
            let mut mass = weights[1];
            let mut residue = increment;
            while residue != 0 {
                budget.charge(1)?;
                table[1][residue] = mass;
                mass = add(mass, weights[1])?;
                residue += increment;
                if residue >= base {
                    residue -= base;
                }
            }
        }
        for i in 2..weights.len() {
            let weight = weights[i];
            let divisor = gcd(weights[0], weight) as usize;
            let (previous, current) = table.split_at_mut(i);
            let previous = &previous[i - 1];
            let current = &mut current[0];
            if weight >= previous[(weight % weights[0]) as usize] {
                budget.charge(base)?;
                current.copy_from_slice(previous);
                continue;
            }
            if divisor == 1 {
                let increment = (weight % weights[0]) as usize;
                let mut value = 0;
                let mut residue = 0;
                budget.charge(base)?;
                for _ in 0..base {
                    value = add(value, weight)?;
                    residue += increment;
                    if residue >= base {
                        residue -= base;
                    }
                    value = value.min(previous[residue]);
                    current[residue] = value;
                }
            } else {
                // Source cache-block traversal. Witness counters are unnecessary
                // for getAllDecompositions and do not affect residue values.
                let mut cur = (weight % weights[0]) as usize;
                let mut prev = 0;
                let increment = cur
                    .checked_sub(divisor)
                    .ok_or_else(|| invalid("invalid residue-table increment"))?;
                budget.charge(divisor)?;
                current[1..divisor].copy_from_slice(&previous[1..divisor]);
                for _ in 1..base / divisor {
                    budget.charge(divisor)?;
                    for _ in 0..divisor {
                        current[cur] = add(current[prev], weight)?.min(previous[cur]);
                        prev += 1;
                        cur += 1;
                    }
                    prev = cur - divisor;
                    cur += increment;
                    if cur >= base {
                        cur -= base;
                    }
                }
                loop {
                    budget.charge(divisor)?;
                    let mut changed = false;
                    prev += 1;
                    cur += 1;
                    for _ in 1..divisor {
                        let value = add(current[prev], weight)?;
                        if value < current[cur] {
                            current[cur] = value;
                            changed = true;
                        }
                        prev += 1;
                        cur += 1;
                    }
                    prev = cur - divisor;
                    cur += increment;
                    if cur >= base {
                        cur -= base;
                    }
                    if !changed {
                        break;
                    }
                }
            }
        }
        Ok(Self {
            weights,
            table,
            lcms,
            mass_in_lcms,
            infinity,
        })
    }
    fn collect(
        &self,
        mass: u64,
        index: usize,
        counts: &mut [u32],
        budget: &mut Budget,
        emit: &mut impl FnMut(&[u32], &mut Budget) -> Result<()>,
    ) -> Result<()> {
        budget.charge(1)?;
        if index == 0 {
            let count = mass / self.weights[0];
            if count * self.weights[0] == mass {
                counts[0] = u32::try_from(count)
                    .map_err(|_| invalid("integer decomposition count overflow"))?;
                emit(counts, budget)?;
            }
            return Ok(());
        }
        let lcm = self.lcms[index];
        let in_lcm = self.mass_in_lcms[index];
        let mut residue = mass % self.weights[0];
        let decrement = self.weights[index] % self.weights[0];
        for i in 0..in_lcm {
            budget.charge(1)?;
            let used = i
                .checked_mul(self.weights[index])
                .ok_or_else(|| invalid("integer mass product overflow"))?;
            if mass < used {
                break;
            }
            counts[index] =
                u32::try_from(i).map_err(|_| invalid("integer decomposition count overflow"))?;
            let minimum = self.table[index - 1][residue as usize];
            if minimum != self.infinity {
                let mut remaining = mass - used;
                while remaining >= minimum {
                    self.collect(remaining, index - 1, counts, budget, emit)?;
                    if remaining < lcm {
                        break;
                    }
                    counts[index] = counts[index]
                        .checked_add(
                            u32::try_from(in_lcm)
                                .map_err(|_| invalid("integer decomposition count overflow"))?,
                        )
                        .ok_or_else(|| invalid("integer decomposition count overflow"))?;
                    remaining -= lcm;
                }
            }
            residue = if residue < decrement {
                residue + self.weights[0] - decrement
            } else {
                residue - decrement
            };
        }
        Ok(())
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}
fn integer(value: f64) -> Result<u64> {
    if !value.is_finite() || !(0.0..=MAX_INTEGER as f64).contains(&value) {
        Err(invalid(
            "scaled integer mass outside exact nonnegative range",
        ))
    } else {
        Ok(value as u64)
    }
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| invalid("integer residue-table overflow"))
}
fn invalid(message: &'static str) -> Error {
    Error::InvalidValue(format!("mass decomposition: {message}"))
}
struct Budget {
    work: usize,
    bytes: usize,
}
impl Default for Budget {
    fn default() -> Self {
        Self {
            work: MAX_DECOMPOSITION_WORK,
            bytes: MAX_DECOMPOSITION_BYTES,
        }
    }
}
impl Budget {
    fn charge(&mut self, amount: usize) -> Result<()> {
        self.work = self
            .work
            .checked_sub(amount)
            .ok_or_else(|| invalid("work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, amount: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(amount)
            .ok_or_else(|| invalid("allocation budget exceeded"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enumerate(solver: &IntegerSolver, mass: u64) -> Vec<Vec<u32>> {
        let mut values = Vec::new();
        solver
            .collect(
                mass,
                solver.weights.len() - 1,
                &mut vec![0; solver.weights.len()],
                &mut Budget::default(),
                &mut |counts, _| {
                    values.push(counts.to_vec());
                    Ok(())
                },
            )
            .unwrap();
        values
    }
    fn brute(weights: &[u64], mass: u64) -> Vec<Vec<u32>> {
        fn visit(weights: &[u64], mass: u64, prefix: &mut Vec<u32>, found: &mut Vec<Vec<u32>>) {
            if weights.is_empty() {
                if mass == 0 {
                    found.push(prefix.clone());
                }
                return;
            }
            for count in 0..=mass / weights[0] {
                prefix.push(count as u32);
                visit(&weights[1..], mass - count * weights[0], prefix, found);
                prefix.pop();
            }
        }
        let mut found = Vec::new();
        visit(weights, mass, &mut Vec::new(), &mut found);
        found
    }
    #[test]
    fn source_residue_table_matches_exhaustive_integer_compositions() {
        for weights in [
            vec![3],
            vec![2, 3],
            vec![6, 9, 20],
            vec![10, 6, 15],
            vec![15, 21, 10, 13],
            vec![7, 4, 9],
            vec![6, 6, 6],
        ] {
            let solver = IntegerSolver::new(weights.clone(), &mut Budget::default()).unwrap();
            for mass in 0..=120 {
                let mut actual = enumerate(&solver, mass);
                let mut expected = brute(&weights, mass);
                actual.sort();
                expected.sort();
                assert_eq!(actual, expected, "weights={weights:?}, mass={mass}");
            }
        }
        let solver = IntegerSolver::new(vec![6, 9, 20], &mut Budget::default()).unwrap();
        assert_eq!(
            enumerate(&solver, 60),
            [
                vec![10, 0, 0],
                vec![7, 2, 0],
                vec![4, 4, 0],
                vec![1, 6, 0],
                vec![0, 0, 3]
            ]
        );
    }

    fn synthetic(masses: &[f64], precision: f64, tolerance: f64) -> MassDecompositionAlgorithm {
        let mut records = Vec::new();
        let mut names = Vec::new();
        for (i, &symbol) in DecompositionResidueSet::Natural19WithoutI
            .symbols()
            .iter()
            .enumerate()
        {
            let name = format!("weight-{symbol}");
            names.push(name.clone());
            records.push(
                ResidueModification::from_record(super::super::ModificationRecord {
                    name,
                    origin: Some(char::from(symbol)),
                    mono_mass: masses.get(i).copied().unwrap_or(1000.0),
                    ..Default::default()
                })
                .unwrap(),
            );
        }
        let registry = ModificationsDB::from_records(records).unwrap();
        MassDecompositionAlgorithm::with_registry(
            MassDecompositionOptions {
                decomp_weights_precision: precision,
                tolerance,
                fixed_modifications: names,
                ..Default::default()
            },
            &registry,
        )
        .unwrap()
    }
    #[test]
    fn source_last_weight_infinity_can_exclude_a_mathematical_composition() {
        // The source uses first*last=12 as its finite "infinity". With these
        // unsorted weights its cache leaves residue1 at12, losing the true13.
        let solver = IntegerSolver::new(vec![6, 9, 2, 2], &mut Budget::default()).unwrap();
        assert_eq!(solver.infinity, 12);
        assert_eq!(
            brute(&solver.weights, 13),
            [vec![0, 1, 0, 2], vec![0, 1, 1, 1], vec![0, 1, 2, 0]]
        );
        assert_eq!(enumerate(&solver, 13), [vec![0, 1, 1, 1], vec![0, 1, 0, 2]]);
    }

    #[test]
    fn rounding_gcd_exclusive_upper_bound_and_zero_tolerance_are_literal() {
        let algorithm = synthetic(&[3.0, 5.0, 8.0], 0.1, 0.0);
        assert_eq!(algorithm.precision, 1.0);
        assert_eq!(&algorithm.solver.weights[..3], [3, 5, 8]);
        assert_eq!(algorithm.minimum_error, 0.0);
        assert_eq!(algorithm.maximum_error, 0.0);
        assert!(algorithm.decompositions(8.0).unwrap().is_empty()); // start=end=8.
        let algorithm = synthetic(&[3.0, 5.0, 8.0], 0.1, 1.0);
        // Search7..9 includes masses7,8, but excludes exact upper boundary9.
        let mut values: Vec<_> = algorithm
            .decompositions(8.0)
            .unwrap()
            .iter()
            .map(|v| v.to_text().unwrap())
            .collect();
        values.sort();
        assert_eq!(values, ["A1 C1", "D1"]);
        let rounded = synthetic(&[2.25, 3.75], 0.5, 0.1);
        assert_eq!(&rounded.solver.weights[..2], [5, 8]);
        assert_eq!(rounded.minimum_error, 0.0);
        assert_eq!(rounded.maximum_error, (0.5 * 5.0 - 2.25) / 2.25);
    }
    #[test]
    fn cumulative_work_output_and_existing_payload_limits_are_atomic() {
        assert!(MassDecomposition::from_counts([(b'A', u32::MAX)].into_iter()).is_err());
        let algorithm = synthetic(&[3.0, 5.0, 8.0], 0.1, 1.0);
        let mut output = vec![MassDecomposition::parse("Z1").unwrap()];
        let original = output.clone();
        let pointer = output.as_ptr();
        let mut budget = Budget {
            work: 40,
            bytes: MAX_DECOMPOSITION_BYTES,
        };
        assert!(
            algorithm
                .append_with_budget(&mut output, 80.0, &mut budget)
                .is_err()
        );
        assert_eq!(output, original);
        assert_eq!(output.as_ptr(), pointer);
        let mut budget = Budget {
            work: MAX_DECOMPOSITION_WORK,
            bytes: 10,
        };
        assert!(
            algorithm
                .append_with_budget(&mut output, 8.0, &mut budget)
                .is_err()
        );
        assert_eq!(output, original);
        let mut output = vec![MassDecomposition::new(); MAX_DECOMPOSITION_OUTPUTS + 1];
        assert!(algorithm.append_decompositions(&mut output, 1.0).is_err());
        assert_eq!(output.len(), MAX_DECOMPOSITION_OUTPUTS + 1);
        let text = (b'A'..=b'Z')
            .map(|c| format!("{}1", char::from(c)))
            .collect::<Vec<_>>()
            .join(" ");
        let large = MassDecomposition::parse(&text).unwrap();
        let mut output = vec![large; MAX_DECOMPOSITION_BYTES / (26 * 64) + 1];
        let old = output.len();
        assert!(algorithm.append_decompositions(&mut output, 1.0).is_err());
        assert_eq!(output.len(), old);
    }
}
