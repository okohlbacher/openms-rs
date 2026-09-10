// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Bounded fine-isotope enumeration of independent multinomial configurations.
//! Source conventions: OpenMS4-core 7c029e8 FineIsotopePatternGenerator and
//! IsoSpecWrapper. This native heap enumerator does not emulate IsoSpec layers.

use super::isotopes::{IsotopeDistribution, IsotopePeak};
use super::{Atom, EmpiricalFormula, element};
use crate::kernel::Peak1D;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashSet};
use std::sync::Arc;

pub const MAX_FINE_ATOMS: usize = 1_000_000;
pub const MAX_FINE_PEAKS: usize = 100_000;
pub const MAX_FINE_STATES: usize = 250_000;
pub const MAX_FINE_BYTES: usize = 128 * 1024 * 1024;
pub const MAX_FINE_WORK: usize = 100_000_000;
/// Includes all custom input categories, even rows with zero atoms.
pub const MAX_FINE_CUSTOM_CELLS: usize = 1_000_000;

/// Probability selection without renormalizing the source isotope abundances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FineIsotopeStop {
    /// Stop when the sum of widened f32 output probabilities reaches 1-value.
    /// Value must be in `[0, 1]`; one requests an empty distribution.
    UnexplainedProbability(f64),
    /// Retain configurations at or above this absolute probability.
    AbsoluteThreshold(f64),
    /// Retain configurations at or above this fraction of the global mode.
    RelativeThreshold(f64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct FineIsotopePatternGenerator {
    pub stop: FineIsotopeStop,
}
impl Default for FineIsotopePatternGenerator {
    fn default() -> Self {
        Self {
            stop: FineIsotopeStop::UnexplainedProbability(0.01),
        }
    }
}

/// One configuration with raw binary64 weight and its finite natural logarithm.
/// Probability may underflow to zero while log_probability remains informative.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FineIsotopeConfiguration {
    pub mass: f64,
    pub probability: f64,
    pub log_probability: f64,
}
impl FineIsotopeConfiguration {
    /// Convert to the source Peak1D representation, narrowing only the final
    /// probability. Public fields are checked; probability/log consistency is
    /// not required because underflow makes that relation lossy.
    pub fn to_peak(&self) -> Result<Peak1D> {
        if self.mass < 0.0 || self.probability < 0.0 {
            return Err(invalid(
                "fine configuration mass and probability must be nonnegative",
            ));
        }
        finite(self.mass)?;
        finite(self.probability)?;
        finite(self.log_probability)?;
        let probability = self.probability as f32;
        finite(f64::from(probability))?;
        Ok(Peak1D::new(self.mass, probability))
    }
}

/// Owning ordered fine-isotope stream. A valid prefix may precede one checked
/// error; errors and normal exhaustion permanently fuse the iterator. Dropping
/// it or stopping after an item does not expand that item's neighbors.
///
/// Streams have no materialized-output cap or full-support preflight. Atom,
/// visited-state, retained-memory and work ceilings bound their entire lifetime.
pub struct FineIsotopeIterator {
    search: Option<ConfigurationSearch>,
    work: FineIsotopeWork,
    cutoff: Cutoff,
    mode_log_probability: f64,
}
impl FineIsotopeIterator {
    /// Enumerate the formula's atom map, ignoring its charge (positive or
    /// negative), as the raw source wrapper does. Natural weights retain source
    /// f32 rounding; fixed labels have unit probability.
    pub fn from_formula(formula: &EmpiricalFormula) -> Result<Self> {
        let mut work = FineIsotopeWork::default();
        validate_formula(formula, 0, &mut work)?;
        let model = Model::from_formula(formula, 0, &mut work)?;
        Self::from_model(model, work)
    }

    /// Copy independent custom populations without normalizing or narrowing
    /// their f64 masses/weights. Rows must be nonempty and shapes match, including
    /// zero-count rows; masses are finite/nonnegative and weights finite/positive.
    /// Empty outer input denotes the single empty configuration.
    pub fn from_isotopes(
        atom_counts: &[u32],
        isotope_masses: &[Vec<f64>],
        isotope_probabilities: &[Vec<f64>],
    ) -> Result<Self> {
        let mut work = FineIsotopeWork::default();
        let model = Model::from_isotopes(
            atom_counts,
            isotope_masses,
            isotope_probabilities,
            &mut work,
        )?;
        Self::from_model(model, work)
    }

    fn from_model(model: Model, mut work: FineIsotopeWork) -> Result<Self> {
        let search = ConfigurationSearch::new(model, &mut work)?;
        Ok(Self {
            mode_log_probability: search.mode_log_probability,
            search: Some(search),
            work,
            cutoff: Cutoff::Log(f64::NEG_INFINITY),
        })
    }

    /// Replace the cutoff for remaining states, comparing raw returned f64
    /// probabilities inclusively. Zero permits even underflowed probabilities.
    /// Unlike the materializer's log comparison, this retains getter equality.
    /// Exhausted streams never restart.
    pub fn with_absolute_threshold(mut self, threshold: f64) -> Result<Self> {
        threshold_log(threshold, None)?;
        self.cutoff = Cutoff::Absolute(threshold);
        Ok(self)
    }

    /// Replace the remaining-state cutoff relative to the ORIGINAL global mode,
    /// even after items have been consumed. Uses raw p/mode_p when both are
    /// finite/positive, otherwise the canonical log cutoff preserves tiny tails.
    /// Values above one exhaust immediately; zero admits all remaining states.
    /// Does not revive exhausted streams.
    pub fn with_relative_threshold(mut self, threshold: f64) -> Result<Self> {
        threshold_log(threshold, None)?;
        self.cutoff = Cutoff::Relative {
            value: threshold,
            mode_log: self.mode_log_probability,
            mode_probability: self.mode_log_probability.exp(),
        };
        Ok(self)
    }
}
impl Iterator for FineIsotopeIterator {
    type Item = Result<FineIsotopeConfiguration>;
    fn next(&mut self) -> Option<Self::Item> {
        let search = self.search.as_mut()?;
        match search.next_configuration(&mut self.work, &self.cutoff, None) {
            Ok(Some(configuration)) => Some(Ok(configuration)),
            Ok(None) => {
                self.search = None;
                None
            }
            Err(error) => {
                self.search = None;
                Some(Err(error))
            }
        }
    }
}
impl std::iter::FusedIterator for FineIsotopeIterator {}

enum Cutoff {
    // Existing materializers retain their established log comparison policy.
    Log(f64),
    Absolute(f64),
    Relative {
        value: f64,
        mode_log: f64,
        mode_probability: f64,
    },
}
impl Cutoff {
    fn includes(&self, log_probability: f64) -> bool {
        match *self {
            Self::Log(value) => log_probability >= value,
            Self::Absolute(value) => value == 0.0 || log_probability.exp() >= value,
            Self::Relative {
                value,
                mode_log,
                mode_probability,
            } => {
                if value > 1.0 {
                    return false;
                }
                if value == 0.0 {
                    return true;
                }
                let probability = log_probability.exp();
                if probability.is_finite()
                    && probability > 0.0
                    && mode_probability.is_finite()
                    && mode_probability > 0.0
                {
                    probability / mode_probability >= value
                } else {
                    log_probability >= mode_log + value.ln()
                }
            }
        }
    }
}

/// Cumulative work allowance shared by theoretical fragment generation.
pub(crate) struct FineIsotopeWork {
    remaining: usize,
}
impl Default for FineIsotopeWork {
    fn default() -> Self {
        Self {
            remaining: MAX_FINE_WORK,
        }
    }
}
impl FineIsotopeWork {
    pub(crate) fn consume(&mut self, units: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(units)
            .ok_or_else(|| invalid("fine isotope work limit exceeded"))?;
        Ok(())
    }
}

#[derive(Clone)]
struct State {
    log_probability: f64,
    counts: Arc<[u32]>,
}
impl PartialEq for State {
    fn eq(&self, other: &Self) -> bool {
        self.log_probability == other.log_probability && self.counts == other.counts
    }
}
impl Eq for State {}
impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        self.log_probability
            .partial_cmp(&other.log_probability)
            .expect("finite configuration log probability")
            .then_with(|| self.counts.cmp(&other.counts))
    }
}

struct Group {
    start: usize,
    end: usize,
    count: u32,
}
#[derive(Default)]
struct Model {
    groups: Vec<Group>,
    masses: Vec<f64>,
    probabilities: Vec<f64>,
    fixed_mass: f64,
    fixed_log_probability: f64,
    custom: bool,
}
impl Model {
    fn from_formula(
        formula: &EmpiricalFormula,
        hydrogen: u32,
        work: &mut FineIsotopeWork,
    ) -> Result<Self> {
        let mut model = Self::default();
        let mut hydrogen_seen = false;
        for (&atom, &count) in &formula.atoms {
            let natural_hydrogen = atom.symbol == "H" && atom.isotope.is_none();
            let count = count as u32
                + if natural_hydrogen {
                    hydrogen_seen = true;
                    hydrogen
                } else {
                    0
                };
            work.consume(element(atom.symbol).unwrap().isotopes().len())?;
            model.add(atom, count)?;
        }
        if !hydrogen_seen && hydrogen > 0 {
            model.add(Atom::resolve("H", None)?, hydrogen)?;
        }
        finite(model.fixed_mass)?;
        finite(model.fixed_log_probability)?;
        Ok(model)
    }

    fn from_isotopes(
        atom_counts: &[u32],
        isotope_masses: &[Vec<f64>],
        isotope_probabilities: &[Vec<f64>],
        work: &mut FineIsotopeWork,
    ) -> Result<Self> {
        if atom_counts.len() != isotope_masses.len()
            || atom_counts.len() != isotope_probabilities.len()
        {
            return Err(invalid("custom isotope outer vector lengths must match"));
        }
        if atom_counts.len() > MAX_FINE_CUSTOM_CELLS {
            return Err(invalid("custom isotope input cell limit exceeded"));
        }
        work.consume(atom_counts.len().saturating_add(1))?;
        let mut cells = 0_usize;
        let mut atoms = 0_usize;
        // Inspect sizes and charge every input cell before walking values or
        // copying anything, even when all population counts are zero.
        for ((&count, masses), probabilities) in atom_counts
            .iter()
            .zip(isotope_masses)
            .zip(isotope_probabilities)
        {
            if masses.is_empty() || masses.len() != probabilities.len() {
                return Err(invalid(
                    "custom isotope rows must be nonempty with matching lengths",
                ));
            }
            cells = cells
                .checked_add(masses.len())
                .filter(|&n| n <= MAX_FINE_CUSTOM_CELLS)
                .ok_or_else(|| invalid("custom isotope input cell limit exceeded"))?;
            atoms = atoms
                .checked_add(count as usize)
                .filter(|&n| n <= MAX_FINE_ATOMS)
                .ok_or_else(|| invalid("fine isotope atom limit exceeded"))?;
        }
        check_memory(
            cells
                .saturating_mul(64)
                .saturating_add(atom_counts.len().saturating_mul(64)),
            0,
            1,
            0,
        )?;
        work.consume(cells.saturating_mul(4))?;
        for (masses, probabilities) in isotope_masses.iter().zip(isotope_probabilities) {
            for (&mass, &probability) in masses.iter().zip(probabilities) {
                if !mass.is_finite() || mass < 0.0 || !probability.is_finite() || probability <= 0.0
                {
                    return Err(invalid(
                        "custom isotope masses must be finite/nonnegative and weights finite/positive",
                    ));
                }
            }
        }
        let mut model = Self {
            custom: true,
            ..Self::default()
        };
        for ((&count, masses), probabilities) in atom_counts
            .iter()
            .zip(isotope_masses)
            .zip(isotope_probabilities)
        {
            if count == 0 {
                continue;
            }
            if masses.len() == 1 {
                model.fixed_mass += f64::from(count) * masses[0];
                model.fixed_log_probability += f64::from(count) * probabilities[0].ln();
            } else {
                let start = model.masses.len();
                reserve(&mut model.masses, masses.len())?;
                reserve(&mut model.probabilities, probabilities.len())?;
                model.masses.extend_from_slice(masses);
                model.probabilities.extend_from_slice(probabilities);
                model.groups.push(Group {
                    start,
                    end: model.masses.len(),
                    count,
                });
            }
        }
        finite(model.fixed_mass)?;
        finite(model.fixed_log_probability)?;
        Ok(model)
    }

    fn add(&mut self, atom: Atom, count: u32) -> Result<()> {
        if count == 0 {
            return Ok(());
        }
        let isotopes = element(atom.symbol)
            .expect("validated formula element")
            .isotopes();
        if let Some(number) = atom.isotope {
            let mass = isotopes
                .iter()
                .find(|i| i.mass_number == number)
                .expect("validated isotope label")
                .mass;
            self.fixed_mass += f64::from(count) * mass;
            return Ok(());
        }
        let start = self.masses.len();
        for isotope in isotopes {
            let probability = f64::from(isotope.abundance as f32);
            if probability > 0.0 {
                self.masses.push(isotope.mass);
                self.probabilities.push(probability);
            }
        }
        let end = self.masses.len();
        if end == start {
            return Err(invalid("natural element has no positive isotope abundance"));
        }
        if end - start == 1 {
            self.fixed_mass += f64::from(count) * self.masses.pop().unwrap();
            self.fixed_log_probability += f64::from(count) * self.probabilities.pop().unwrap().ln();
        } else {
            self.groups.push(Group { start, end, count });
        }
        Ok(())
    }
}

impl FineIsotopePatternGenerator {
    /// Return distinct configurations sorted by mass, with source f32 output
    /// rounding and no renormalization. Positive charge adds natural H atoms;
    /// negative charge errors. Masses are not divided by charge.
    ///
    /// Inclusive thresholds and deterministic configuration ties are native
    /// policies. Coverage accumulates widened f32 probabilities, not raw weights.
    /// Full zero-threshold support includes output-underflow zeros. Limits return
    /// errors without a partial distribution; see module constants.
    pub fn run(&self, formula: &EmpiricalFormula) -> Result<IsotopeDistribution> {
        self.run_with_work(formula, &mut FineIsotopeWork::default())
    }

    /// Materialize custom independent isotope populations using this generator's
    /// selection policy. Input f64 weights are preserved without normalization;
    /// only final output probabilities are narrowed to f32. Validation also
    /// applies to zero-count rows and requests for zero coverage.
    pub fn run_with_isotopes(
        &self,
        atom_counts: &[u32],
        isotope_masses: &[Vec<f64>],
        isotope_probabilities: &[Vec<f64>],
    ) -> Result<IsotopeDistribution> {
        let selection = Selection::new(self.stop)?;
        let mut work = FineIsotopeWork::default();
        let model = Model::from_isotopes(
            atom_counts,
            isotope_masses,
            isotope_probabilities,
            &mut work,
        )?;
        materialize(model, selection, &mut work)
    }

    pub(crate) fn run_with_work(
        &self,
        formula: &EmpiricalFormula,
        work: &mut FineIsotopeWork,
    ) -> Result<IsotopeDistribution> {
        let selection = Selection::new(self.stop)?;
        if formula.charge < 0 {
            return Err(invalid(
                "fine isotopes do not support negative formula charge",
            ));
        }
        validate_formula(formula, formula.charge as u32, work)?;
        if selection.coverage == Some(0.0) {
            return Ok(IsotopeDistribution::empty());
        }
        let model = Model::from_formula(formula, formula.charge as u32, work)?;
        materialize(model, selection, work)
    }
}

struct Selection {
    coverage: Option<f64>,
    threshold: f64,
    relative: bool,
}
impl Selection {
    fn new(stop: FineIsotopeStop) -> Result<Self> {
        match stop {
            FineIsotopeStop::UnexplainedProbability(value) => {
                if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                    return Err(invalid("unexplained probability must lie in [0,1]"));
                }
                Ok(Self {
                    coverage: Some(1.0 - value),
                    threshold: 0.0,
                    relative: false,
                })
            }
            FineIsotopeStop::AbsoluteThreshold(value)
            | FineIsotopeStop::RelativeThreshold(value) => {
                threshold_log(value, None)?;
                Ok(Self {
                    coverage: None,
                    threshold: value,
                    relative: matches!(stop, FineIsotopeStop::RelativeThreshold(_)),
                })
            }
        }
    }
}

fn materialize(
    model: Model,
    selection: Selection,
    work: &mut FineIsotopeWork,
) -> Result<IsotopeDistribution> {
    if selection.coverage == Some(0.0) {
        return Ok(IsotopeDistribution::empty());
    }
    // This preflight belongs only to materialization. A stream may consume a
    // small prefix even when the full support is too large to collect.
    if selection.coverage.is_none() && selection.threshold == 0.0 {
        let mut total = 1_usize;
        for group in &model.groups {
            total = total
                .checked_mul(support_count(group.count as usize, group.end - group.start))
                .filter(|&n| n <= MAX_FINE_PEAKS)
                .ok_or_else(|| invalid("full fine isotope support exceeds peak limit"))?;
        }
    }
    let mut search = ConfigurationSearch::new(model, work)?;
    let cutoff = Cutoff::Log(threshold_log(
        selection.threshold,
        selection.relative.then_some(search.mode_log_probability),
    )?);
    if !cutoff.includes(search.mode_log_probability) {
        return Ok(IsotopeDistribution::empty());
    }
    let mut output = Vec::new();
    let mut accumulated = 0.0;
    while let Some(configuration) = search.next_configuration(work, &cutoff, Some(output.len()))? {
        let peak = configuration.to_peak()?;
        let probability = f64::from(peak.intensity);
        reserve(&mut output, 1)?;
        output.push(IsotopePeak {
            mass: peak.mz,
            probability,
        });
        accumulated = finite(accumulated + probability)?;
        if selection
            .coverage
            .is_some_and(|target| accumulated >= target)
        {
            break;
        }
    }
    work.consume(heap_work(output.len(), 1).saturating_mul(output.len()))?;
    output.sort_by(|a, b| a.mass.total_cmp(&b.mass));
    IsotopeDistribution::from_peaks(output)
}

struct ConfigurationSearch {
    model: Model,
    factorials: Vec<f64>,
    logs: Vec<f64>,
    factorial_constant: f64,
    mode_log_probability: f64,
    base_bytes: usize,
    seen: HashSet<Arc<[u32]>>,
    heap: BinaryHeap<State>,
    pending: Option<State>,
}
impl ConfigurationSearch {
    fn new(model: Model, work: &mut FineIsotopeWork) -> Result<Self> {
        finite(model.fixed_mass)?;
        finite(model.fixed_log_probability)?;
        let dimension = model.masses.len();
        let maximum_count = model
            .groups
            .iter()
            .map(|g| g.count as usize)
            .max()
            .unwrap_or(0);
        let base_bytes = (maximum_count + 1)
            .saturating_mul(8)
            .saturating_add(dimension.saturating_mul(64))
            .saturating_add(model.groups.len().saturating_mul(64));
        check_memory(base_bytes, dimension, 1, 0)?;
        work.consume(
            maximum_count
                .saturating_add(dimension.saturating_mul(4))
                .saturating_add(1),
        )?;
        let mut factorials = Vec::new();
        reserve(&mut factorials, maximum_count + 1)?;
        factorials.push(0.0);
        for value in 1..=maximum_count {
            factorials.push(factorials[value - 1] + (value as f64).ln());
        }
        let logs: Vec<_> = model.probabilities.iter().map(|p| p.ln()).collect();
        let factorial_constant = model.fixed_log_probability
            + compensated(model.groups.iter().map(|g| factorials[g.count as usize]));
        let mut mode = vec![0; dimension];
        for group in &model.groups {
            let probabilities = &model.probabilities[group.start..group.end];
            if model.custom {
                // Scale only mode-finding weights. Keep raw f64 values and logs
                // intact; very small scaled values may underflow to zero, but
                // their original categories still participate in enumeration.
                work.consume(probabilities.len().saturating_mul(3))?;
                let maximum = probabilities.iter().copied().fold(0.0, f64::max);
                let scaled: Vec<_> = probabilities.iter().map(|p| p / maximum).collect();
                find_mode(
                    group.count,
                    &scaled,
                    &mut mode[group.start..group.end],
                    work,
                )?;
            } else {
                find_mode(
                    group.count,
                    probabilities,
                    &mut mode[group.start..group.end],
                    work,
                )?;
            }
        }
        let mode_log_probability = log_probability(&mode, &logs, &factorials, factorial_constant)?;
        let counts: Arc<[u32]> = mode.into();
        let mut seen = HashSet::new();
        seen.try_reserve(1)
            .map_err(|_| invalid("fine isotope allocation failed"))?;
        seen.insert(Arc::clone(&counts));
        let mut heap = BinaryHeap::new();
        heap.try_reserve(1)
            .map_err(|_| invalid("fine isotope allocation failed"))?;
        heap.push(State {
            log_probability: mode_log_probability,
            counts,
        });
        Ok(Self {
            model,
            factorials,
            logs,
            factorial_constant,
            mode_log_probability,
            base_bytes,
            seen,
            heap,
            pending: None,
        })
    }

    fn next_configuration(
        &mut self,
        work: &mut FineIsotopeWork,
        cutoff: &Cutoff,
        retained_peaks: Option<usize>,
    ) -> Result<Option<FineIsotopeConfiguration>> {
        if !cutoff.includes(self.mode_log_probability)
            || self
                .pending
                .as_ref()
                .is_some_and(|state| !cutoff.includes(state.log_probability))
        {
            return Ok(None);
        }
        if let Some(state) = self.pending.take() {
            self.expand(&state, work, retained_peaks.unwrap_or(0))?;
        }
        if self.heap.is_empty() {
            return Ok(None);
        }
        let dimension = self.model.masses.len();
        work.consume(heap_work(self.heap.len(), dimension))?;
        let state = self.heap.pop().expect("nonempty configuration frontier");
        if !cutoff.includes(state.log_probability) {
            return Ok(None);
        }
        if retained_peaks.is_some_and(|n| n >= MAX_FINE_PEAKS) {
            return Err(invalid("fine isotope output peak limit exceeded"));
        }
        check_memory(
            self.base_bytes,
            dimension,
            self.seen.len(),
            retained_peaks.map_or(0, |n| n + 1),
        )?;
        work.consume(dimension.saturating_mul(2).saturating_add(1))?;
        let probability = finite(state.log_probability.exp())?;
        let mass = finite(
            self.model.fixed_mass
                + compensated(
                    state
                        .counts
                        .iter()
                        .zip(&self.model.masses)
                        .map(|(&count, &mass)| f64::from(count) * mass),
                ),
        )?;
        let configuration = FineIsotopeConfiguration {
            mass,
            probability,
            log_probability: state.log_probability,
        };
        // Coverage can stop here without paying for or failing in expansion of
        // the final returned configuration. Expansion resumes on the next call.
        self.pending = Some(state);
        Ok(Some(configuration))
    }

    fn expand(
        &mut self,
        state: &State,
        work: &mut FineIsotopeWork,
        retained_peaks: usize,
    ) -> Result<()> {
        let dimension = self.model.masses.len();
        work.consume(dimension)?;
        let mut candidate = state.counts.to_vec();
        // Include all transfers and plateaus. Multinomial superlevel sets have
        // nonincreasing paths from a mode to every configuration.
        for group in &self.model.groups {
            for donor in group.start..group.end {
                if candidate[donor] == 0 {
                    continue;
                }
                for receiver in group.start..group.end {
                    if donor == receiver {
                        continue;
                    }
                    work.consume(dimension.saturating_add(1))?;
                    candidate[donor] -= 1;
                    candidate[receiver] += 1;
                    if !self.seen.contains(candidate.as_slice()) {
                        if self.seen.len() >= MAX_FINE_STATES {
                            return Err(invalid("fine isotope visited-state limit exceeded"));
                        }
                        check_memory(
                            self.base_bytes,
                            dimension,
                            self.seen.len() + 1,
                            retained_peaks,
                        )?;
                        work.consume(
                            dimension
                                .saturating_mul(4)
                                .saturating_add(heap_work(self.heap.len() + 1, dimension)),
                        )?;
                        let log_probability = log_probability(
                            &candidate,
                            &self.logs,
                            &self.factorials,
                            self.factorial_constant,
                        )?;
                        let counts: Arc<[u32]> = candidate.clone().into();
                        self.seen
                            .try_reserve(1)
                            .map_err(|_| invalid("fine isotope allocation failed"))?;
                        self.heap
                            .try_reserve(1)
                            .map_err(|_| invalid("fine isotope allocation failed"))?;
                        self.seen.insert(Arc::clone(&counts));
                        self.heap.push(State {
                            log_probability,
                            counts,
                        });
                    }
                    candidate[donor] += 1;
                    candidate[receiver] -= 1;
                }
            }
        }
        Ok(())
    }
}

fn validate_formula(
    formula: &EmpiricalFormula,
    hydrogen: u32,
    work: &mut FineIsotopeWork,
) -> Result<()> {
    work.consume(formula.atoms.len().saturating_add(1))?;
    let mut total = hydrogen as usize;
    for &count in formula.atoms.values() {
        if count < 0 {
            return Err(invalid("fine isotopes require nonnegative atom counts"));
        }
        total = total
            .checked_add(count as usize)
            .filter(|&n| n <= MAX_FINE_ATOMS)
            .ok_or_else(|| invalid("fine isotope atom limit exceeded"))?;
    }
    if total > MAX_FINE_ATOMS {
        return Err(invalid("fine isotope atom limit exceeded"));
    }
    Ok(())
}

fn threshold_log(value: f64, mode: Option<f64>) -> Result<f64> {
    if !value.is_finite() || value < 0.0 {
        return Err(invalid(
            "fine isotope threshold must be finite and nonnegative",
        ));
    }
    Ok(if value == 0.0 {
        f64::NEG_INFINITY
    } else {
        mode.unwrap_or(0.0) + value.ln()
    })
}

fn find_mode(
    n: u32,
    probabilities: &[f64],
    counts: &mut [u32],
    work: &mut FineIsotopeWork,
) -> Result<()> {
    let sum: f64 = probabilities.iter().sum();
    let mut allocated = 0_u32;
    for (count, &probability) in counts.iter_mut().zip(probabilities) {
        *count = (f64::from(n) * probability / sum).floor() as u32;
        allocated += *count;
    }
    // Floors select a prefix of globally largest marginal gains. Greedy
    // allocation of the (normally < isotope count) remainder gives a mode.
    while allocated != n {
        work.consume(counts.len())?;
        let mut best = None;
        for i in 0..counts.len() {
            if allocated > n && counts[i] == 0 {
                continue;
            }
            if best.is_none_or(|j| {
                if allocated < n {
                    probabilities[i] * f64::from(counts[j] + 1)
                        > probabilities[j] * f64::from(counts[i] + 1)
                } else {
                    probabilities[i] * f64::from(counts[j])
                        < probabilities[j] * f64::from(counts[i])
                }
            }) {
                best = Some(i);
            }
        }
        let best = best.expect("positive isotope support");
        if allocated < n {
            counts[best] += 1;
            allocated += 1;
        } else {
            counts[best] -= 1;
            allocated -= 1;
        }
    }
    // Verify the local exchange certificate after floating-point initialization.
    // With bounded counts and f32 abundances the comparison products fit f64.
    loop {
        let mut transfer = None;
        'search: for i in 0..counts.len() {
            if counts[i] == 0 {
                continue;
            }
            for j in 0..counts.len() {
                work.consume(1)?;
                if i != j
                    && f64::from(counts[i]) * probabilities[j]
                        > f64::from(counts[j] + 1) * probabilities[i]
                {
                    transfer = Some((i, j));
                    break 'search;
                }
            }
        }
        let Some((i, j)) = transfer else {
            break;
        };
        counts[i] -= 1;
        counts[j] += 1;
    }
    Ok(())
}
fn log_probability(counts: &[u32], logs: &[f64], factorials: &[f64], constant: f64) -> Result<f64> {
    finite(
        constant
            + compensated(
                counts
                    .iter()
                    .zip(logs)
                    .map(|(&n, &p)| f64::from(n) * p - factorials[n as usize]),
            ),
    )
}
fn compensated(values: impl Iterator<Item = f64>) -> f64 {
    let mut sum: f64 = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}
fn support_count(n: usize, isotopes: usize) -> usize {
    let mut value = 1_usize;
    for i in 1..isotopes {
        let Some(product) = value.checked_mul(n + i) else {
            return MAX_FINE_PEAKS + 1;
        };
        value = product / i;
        if value > MAX_FINE_PEAKS {
            return MAX_FINE_PEAKS + 1;
        }
    }
    value
}
fn heap_work(length: usize, dimension: usize) -> usize {
    ((usize::BITS - length.leading_zeros()) as usize + 1).saturating_mul(dimension + 1)
}
fn check_memory(base: usize, dimension: usize, states: usize, peaks: usize) -> Result<()> {
    // Covers shared state payloads, both container capacities/headers, candidate
    // scratch, factorials and geometric output capacity; not merely peak count.
    let bytes = base
        .saturating_add(states.saturating_mul(dimension.saturating_mul(8).saturating_add(192)))
        .saturating_add(peaks.saturating_mul(64));
    if bytes > MAX_FINE_BYTES {
        Err(invalid("fine isotope retained-payload limit exceeded"))
    } else {
        Ok(())
    }
}
fn reserve<T>(values: &mut Vec<T>, additional: usize) -> Result<()> {
    values
        .try_reserve(additional)
        .map_err(|_| invalid("fine isotope allocation failed"))
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("fine isotope numerical overflow"))
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_weight_modes_keep_both_members_of_a_probability_plateau() {
        let mut counts = [0; 2];
        find_mode(3, &[0.5, 0.5], &mut counts, &mut FineIsotopeWork::default()).unwrap();
        assert_eq!(counts, [2, 1]);
        let factorials = [0.0, 0.0, 2.0_f64.ln(), 6.0_f64.ln()];
        let logs = [0.5_f64.ln(); 2];
        let left = log_probability(&[2, 1], &logs, &factorials, factorials[3]).unwrap();
        let right = log_probability(&[1, 2], &logs, &factorials, factorials[3]).unwrap();
        assert_eq!(left, right);
        assert!((left.exp() - 0.375).abs() < 1e-15);
        let mut seen = HashSet::new();
        seen.insert(Arc::<[u32]>::from([2, 1]));
        assert!(seen.insert(Arc::from([1, 2])));
        // Probability equality does not collapse configuration identity. The
        // heap tie-breaker is deterministic while retaining both tied states.
        let mut heap = BinaryHeap::from([
            State {
                log_probability: left,
                counts: Arc::from([2, 1]),
            },
            State {
                log_probability: right,
                counts: Arc::from([1, 2]),
            },
        ]);
        assert_eq!(&*heap.pop().unwrap().counts, &[2, 1]);
        assert_eq!(&*heap.pop().unwrap().counts, &[1, 2]);
    }

    #[test]
    fn shared_work_is_consumed_across_calls() {
        let formula = EmpiricalFormula::parse("C6H12O6").unwrap();
        let generator = FineIsotopePatternGenerator::default();
        let mut work = FineIsotopeWork::default();
        generator.run_with_work(&formula, &mut work).unwrap();
        let used = MAX_FINE_WORK - work.remaining;
        assert!(used > 0);
        work.remaining = used - 1;
        assert!(generator.run_with_work(&formula, &mut work).is_err());
    }
}
