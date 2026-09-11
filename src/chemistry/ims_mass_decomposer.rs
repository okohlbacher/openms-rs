// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Checked public equivalents of the source IMS integer/real decomposers.
//! Supplied weight order, residue-table sentinels, witness choices and the real
//! wrapper's exclusive upper endpoint are retained. See IMS_DECOMPOSER_SUPPORT.md.

use super::IMSWeights;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::mem::size_of;

pub const MAX_IMS_DECOMPOSER_DIMENSIONS: usize = 128;
pub const MAX_IMS_DECOMPOSER_TABLE_CELLS: usize = 4_000_000;
pub const MAX_IMS_DECOMPOSER_OUTPUTS: usize = 100_000;
pub const MAX_IMS_DECOMPOSER_WORK: usize = 50_000_000;
pub const MAX_IMS_DECOMPOSER_BYTES: usize = 64 * 1024 * 1024;
const U64_END: f64 = 18_446_744_073_709_551_616.0;

/// Native counterpart of the source abstract MassDecomposer interface.
/// Integers use portable u64 masses/u32 multiplicities. Absence of a single
/// decomposition is `None`, distinct from a present all-zero decomposition.
pub trait IMSMassDecomposer {
    fn exists(&self, mass: u64) -> Result<bool>;
    fn decomposition(&self, mass: u64) -> Result<Option<Vec<u32>>>;
    fn decompositions(&self, mass: u64) -> Result<Vec<Vec<u32>>>;
    fn number_of_decompositions(&self, mass: u64) -> Result<u32>;
}

/// An owned, immutable extended residue table. Construction neither sorts nor
/// divides weights by their GCD. Queries have fresh cumulative budgets.
#[derive(Clone, Debug)]
pub struct IMSIntegerMassDecomposer {
    weights: Vec<u64>,
    table: Vec<Vec<u64>>,
    witness: Vec<(usize, u32)>,
    lcms: Vec<u64>,
    mass_in_lcms: Vec<u64>,
    infinity: u64,
}

impl IMSIntegerMassDecomposer {
    pub fn new(weights: &IMSWeights) -> Result<Self> {
        Self::from_integer_weights(weights.weights())
    }
    /// Native convenience for exact integer inputs without an f64 conversion.
    pub fn from_integer_weights(weights: &[u64]) -> Result<Self> {
        Self::build(weights, &mut Work::default())
    }
    fn build(input: &[u64], work: &mut Work) -> Result<Self> {
        let n = input.len();
        if n == 0 || n > MAX_IMS_DECOMPOSER_DIMENSIONS {
            return Err(invalid("alphabet must contain 1..=128 weights"));
        }
        work.consume(n)?;
        if input.contains(&0) {
            return Err(invalid("integer weights must be positive"));
        }
        let mut weights = vector(n, work)?;
        weights.extend_from_slice(input);
        // Source queries other than enumeration have no table for a singleton.
        // The native singleton extension needs neither a table nor a sentinel.
        if n == 1 {
            return Ok(Self {
                weights,
                table: Vec::new(),
                witness: Vec::new(),
                lcms: Vec::new(),
                mass_in_lcms: Vec::new(),
                infinity: 0,
            });
        }
        let base = usize::try_from(input[0]).map_err(|_| invalid("table dimension overflow"))?;
        let cells = base
            .checked_mul(n)
            .filter(|&c| c <= MAX_IMS_DECOMPOSER_TABLE_CELLS)
            .ok_or_else(|| invalid("residue table cell limit exceeded"))?;
        work.consume(cells)?;
        let infinity = product(input[0], input[n - 1])?;
        let mut table = vector(n, work)?;
        for _ in 0..n {
            let mut row = filled(base, infinity, work)?;
            row[0] = 0;
            table.push(row);
        }
        let mut witness = filled(base, (0, 0), work)?;
        let mut lcms = filled(n, 0, work)?;
        let mut mass_in_lcms = filled(n, 0, work)?;
        let increment = (input[1] % input[0]) as usize;
        let mut mass = input[1];
        let mut residue = increment;
        let mut counter = 0u32;
        while residue != 0 {
            work.consume(1)?;
            table[1][residue] = mass;
            mass = add(mass, input[1])?;
            counter = increment_count(counter)?;
            witness[residue] = (1, counter);
            residue += increment;
            if residue >= base {
                residue -= base;
            }
        }
        let divisor = gcd(input[0], input[1], work)?;
        lcms[1] = product(input[1], input[0])? / divisor;
        mass_in_lcms[1] = input[0] / divisor;
        for i in 2..n {
            let weight = input[i];
            let divisor = gcd(input[0], weight, work)? as usize;
            lcms[i] = product(weight, input[0])? / divisor as u64;
            mass_in_lcms[i] = input[0] / divisor as u64;
            let (previous, current) = table.split_at_mut(i);
            let previous = &previous[i - 1];
            let current = &mut current[0];
            if weight >= previous[(weight % input[0]) as usize] {
                work.consume(base)?;
                current.copy_from_slice(previous);
                continue;
            }
            if divisor == 1 {
                let increment = (weight % input[0]) as usize;
                let mut value = 0;
                let mut residue = 0;
                let mut counter = 0u32;
                work.consume(base)?;
                for _ in 0..base {
                    value = add(value, weight)?;
                    residue += increment;
                    counter = increment_count(counter)?;
                    if residue >= base {
                        residue -= base;
                    }
                    if value > previous[residue] {
                        value = previous[residue];
                        counter = 0;
                    } else {
                        witness[residue] = (i, counter);
                    }
                    current[residue] = value;
                }
            } else {
                // Preserve source cache traversal and witness counters, including
                // the single counter increment before its second inner loop.
                let mut cur = (weight % input[0]) as usize;
                let mut prev = 0;
                let increment = cur
                    .checked_sub(divisor)
                    .ok_or_else(|| invalid("invalid residue-table increment"))?;
                let mut counters = filled(base, 0u32, work)?;
                work.consume(divisor)?;
                current[1..divisor].copy_from_slice(&previous[1..divisor]);
                for _ in 1..base / divisor {
                    work.consume(divisor)?;
                    for _ in 0..divisor {
                        counters[cur] = increment_count(counters[cur])?;
                        let value = add(current[prev], weight)?;
                        if value > previous[cur] {
                            current[cur] = previous[cur];
                            counters[cur] = 0;
                        } else {
                            current[cur] = value;
                            witness[cur] = (i, counters[cur]);
                        }
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
                    work.consume(divisor)?;
                    let mut changed = false;
                    prev += 1;
                    cur += 1;
                    counters[cur] = increment_count(counters[cur])?;
                    for _ in 1..divisor {
                        let value = add(current[prev], weight)?;
                        if value < current[cur] {
                            current[cur] = value;
                            changed = true;
                            witness[cur] = (i, counters[cur]);
                        } else {
                            counters[cur] = 0;
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
            witness,
            lcms,
            mass_in_lcms,
            infinity,
        })
    }

    /// This is the source final-table test, which can disagree with enumeration
    /// for certain unsorted alphabets because its infinity sentinel is finite.
    pub fn exists(&self, mass: u64) -> Result<bool> {
        if self.weights.len() == 1 {
            return Ok(mass % self.weights[0] == 0);
        }
        let value = self.table[self.weights.len() - 1][(mass % self.weights[0]) as usize];
        Ok(value != self.infinity && mass >= value)
    }
    /// Follow the source witness vector, rather than selecting the first entry
    /// from enumeration. Invalid/nonprogressing witnesses return a checked error.
    pub fn decomposition(&self, mass: u64) -> Result<Option<Vec<u32>>> {
        let mut work = Work::default();
        self.one_with_work(mass, &mut work)
    }
    fn one_with_work(&self, mass: u64, work: &mut Work) -> Result<Option<Vec<u32>>> {
        if !self.exists(mass)? {
            return Ok(None);
        }
        let mut result = filled(self.weights.len(), 0u32, work)?;
        if self.weights.len() == 1 {
            result[0] = count(mass / self.weights[0])?;
            return Ok(Some(result));
        }
        let mut residue = (mass % self.weights[0]) as usize;
        let mut remaining = self.table[self.weights.len() - 1][residue];
        result[0] = count((mass - remaining) / self.weights[0])?;
        while remaining != 0 {
            work.consume(1)?;
            let (index, number) = self.witness[residue];
            result[index] = result[index]
                .checked_add(number)
                .ok_or_else(|| invalid("decomposition multiplicity overflow"))?;
            let used = product(u64::from(number), self.weights[index])?;
            // Preserve the source's finite early-break branch after adding the
            // witness. Do not silently replace its result with a different one.
            if remaining < used {
                break;
            }
            if used == 0 {
                return Err(invalid("nonprogressing decomposition witness"));
            }
            remaining -= used;
            residue = (remaining % self.weights[0]) as usize;
        }
        Ok(Some(result))
    }
    pub fn decompositions(&self, mass: u64) -> Result<Vec<Vec<u32>>> {
        let mut work = Work::default();
        self.all_with_work(mass, &mut work)
    }
    fn all_with_work(&self, mass: u64, work: &mut Work) -> Result<Vec<Vec<u32>>> {
        let mut output = Vec::new();
        let mut counts = filled(self.weights.len(), 0u32, work)?;
        self.walk(
            mass,
            self.weights.len() - 1,
            &mut counts,
            work,
            &mut |counts, work| retain(&mut output, counts, work),
        )?;
        Ok(output)
    }
    /// Count the same traversal without retaining every decomposition. The
    /// 100,000 materialized-output cap does not apply; work/overflow guards do.
    pub fn number_of_decompositions(&self, mass: u64) -> Result<u32> {
        let mut work = Work::default();
        let mut counts = filled(self.weights.len(), 0u32, &mut work)?;
        let mut total = 0u32;
        self.walk(
            mass,
            self.weights.len() - 1,
            &mut counts,
            &mut work,
            &mut |_, work| {
                work.consume(1)?;
                total = increment_count(total)?;
                Ok(())
            },
        )?;
        Ok(total)
    }
    fn walk(
        &self,
        mass: u64,
        index: usize,
        counts: &mut [u32],
        work: &mut Work,
        emit: &mut impl FnMut(&[u32], &mut Work) -> Result<()>,
    ) -> Result<()> {
        work.consume(1)?;
        if index == 0 {
            let number = mass / self.weights[0];
            if number * self.weights[0] == mass {
                counts[0] = count(number)?;
                emit(counts, work)?;
            }
            return Ok(());
        }
        let lcm = self.lcms[index];
        let in_lcm = self.mass_in_lcms[index];
        let mut residue = mass % self.weights[0];
        let decrement = self.weights[index] % self.weights[0];
        for i in 0..in_lcm {
            work.consume(1)?;
            let used = product(i, self.weights[index])?;
            if mass < used {
                break;
            }
            counts[index] = count(i)?;
            let minimum = self.table[index - 1][residue as usize];
            if minimum != self.infinity {
                let mut remaining = mass - used;
                while remaining >= minimum {
                    self.walk(remaining, index - 1, counts, work, emit)?;
                    if remaining < lcm {
                        break;
                    }
                    counts[index] = counts[index]
                        .checked_add(count(in_lcm)?)
                        .ok_or_else(|| invalid("decomposition multiplicity overflow"))?;
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

impl IMSMassDecomposer for IMSIntegerMassDecomposer {
    fn exists(&self, mass: u64) -> Result<bool> {
        Self::exists(self, mass)
    }
    fn decomposition(&self, mass: u64) -> Result<Option<Vec<u32>>> {
        Self::decomposition(self, mass)
    }
    fn decompositions(&self, mass: u64) -> Result<Vec<Vec<u32>>> {
        Self::decompositions(self, mass)
    }
    fn number_of_decompositions(&self, mass: u64) -> Result<u32> {
        Self::number_of_decompositions(self, mass)
    }
}

/// Real-mass wrapper around the ordered integer traversal. Bounds and filtering
/// use the original masses and supplied precision, including any prior GCD
/// division performed by the caller. There is no automatic GCD transformation.
#[derive(Clone, Debug)]
pub struct IMSRealMassDecomposer {
    weights: IMSWeights,
    integer: IMSIntegerMassDecomposer,
    precision: f64,
    minimum_error: f64,
    maximum_error: f64,
}
impl IMSRealMassDecomposer {
    pub fn new(weights: &IMSWeights) -> Result<Self> {
        let mut work = Work::default();
        if weights.len() > MAX_IMS_DECOMPOSER_DIMENSIONS {
            return Err(invalid("alphabet dimension limit exceeded"));
        }
        // Precharge the two rounding-bound traversals and the owned pair copy.
        work.consume(weights.len() * 6)?;
        work.allocate(weights.len() * 2 * size_of::<u64>())?;
        let minimum_error = weights.min_rounding_error()?;
        let maximum_error = weights.max_rounding_error()?;
        let precision = weights
            .precision()
            .ok_or_else(|| invalid("precision is unset"))?;
        if !precision.is_finite() || precision == 0.0 {
            return Err(invalid(
                "real decomposition requires finite nonzero precision",
            ));
        }
        let integer = IMSIntegerMassDecomposer::build(weights.weights(), &mut work)?;
        Ok(Self {
            weights: weights.clone(),
            integer,
            precision,
            minimum_error,
            maximum_error,
        })
    }
    pub fn decompositions(&self, mass: f64, error: f64) -> Result<Vec<Vec<u32>>> {
        self.decompositions_with_constraints(mass, error, &BTreeMap::new())
    }
    /// Inclusive multiplicity constraints, in input index space. Invalid indices
    /// are checked errors; reversed lower/upper bounds naturally exclude results.
    pub fn decompositions_with_constraints(
        &self,
        mass: f64,
        error: f64,
        constraints: &BTreeMap<usize, (u32, u32)>,
    ) -> Result<Vec<Vec<u32>>> {
        let mut work = Work::default();
        self.all_with_work(mass, error, constraints, &mut work)
    }
    fn all_with_work(
        &self,
        mass: f64,
        error: f64,
        constraints: &BTreeMap<usize, (u32, u32)>,
        work: &mut Work,
    ) -> Result<Vec<Vec<u32>>> {
        work.consume(constraints.len())?;
        if constraints.keys().any(|&index| index >= self.weights.len()) {
            return Err(invalid("constraint index out of range"));
        }
        let (start, end) = self.range(mass, error, false, work)?;
        let mut output = Vec::new();
        let mut counts = filled(self.weights.len(), 0u32, work)?;
        for integer_mass in start..end {
            work.consume(1)?;
            self.integer.walk(
                integer_mass,
                self.weights.len() - 1,
                &mut counts,
                work,
                &mut |counts, work| {
                    if self.accepts(counts, mass, error, work)? {
                        work.consume(constraints.len())?;
                        if constraints
                            .iter()
                            .all(|(&index, &(lo, hi))| counts[index] >= lo && counts[index] <= hi)
                        {
                            retain(&mut output, counts, work)?;
                        }
                    }
                    Ok(())
                },
            )?;
        }
        Ok(output)
    }
    /// Preserve the source count wrapper's distinct start=1 when mass-error<=0.
    /// Thus this need not equal `decompositions(...).len()` around zero mass.
    pub fn number_of_decompositions(&self, mass: f64, error: f64) -> Result<u64> {
        let mut work = Work::default();
        let (start, end) = self.range(mass, error, true, &mut work)?;
        let mut counts = filled(self.weights.len(), 0u32, &mut work)?;
        let mut total = 0u64;
        for integer_mass in start..end {
            work.consume(1)?;
            self.integer.walk(
                integer_mass,
                self.weights.len() - 1,
                &mut counts,
                &mut work,
                &mut |counts, work| {
                    if self.accepts(counts, mass, error, work)? {
                        total = total
                            .checked_add(1)
                            .ok_or_else(|| invalid("decomposition total overflow"))?;
                    }
                    Ok(())
                },
            )?;
        }
        Ok(total)
    }
    fn range(&self, mass: f64, error: f64, counting: bool, work: &mut Work) -> Result<(u64, u64)> {
        if !mass.is_finite() || !error.is_finite() {
            return Err(invalid("mass and error must be finite"));
        }
        let start = if counting && mass - error <= 0.0 {
            1
        } else {
            integer(((1.0 + self.minimum_error) * (mass - error) / self.precision).ceil())?
        };
        let end = integer(((1.0 + self.maximum_error) * (mass + error) / self.precision).floor())?;
        // Charge every candidate integer even when all are later pruned. Avoid a
        // huge empty-range loop independently of the number of retained rows.
        let width = usize::try_from(end.saturating_sub(start))
            .map_err(|_| invalid("integer range exceeds work limit"))?;
        work.consume(width)?;
        Ok((start, end))
    }
    fn accepts(&self, counts: &[u32], mass: f64, error: f64, work: &mut Work) -> Result<bool> {
        // The existing value helper has its own bounded loop; charge it here to
        // keep every candidate inside this operation's cumulative budget.
        work.consume(self.weights.len() * 3)?;
        let parent = self.weights.parent_mass(counts)?;
        Ok((parent - mass).abs() <= error)
    }
}

fn retain(output: &mut Vec<Vec<u32>>, counts: &[u32], work: &mut Work) -> Result<()> {
    if output.len() >= MAX_IMS_DECOMPOSER_OUTPUTS {
        return Err(invalid("materialized output limit exceeded"));
    }
    if output.len() == output.capacity() {
        work.consume(output.len())?; // Moving existing row descriptors on growth.
        let capacity = (output.capacity().max(4) * 2).min(MAX_IMS_DECOMPOSER_OUTPUTS);
        work.allocate(capacity * size_of::<Vec<u32>>())?;
        output
            .try_reserve_exact(capacity - output.len())
            .map_err(|_| invalid("output allocation failed"))?;
    }
    let mut row = vector(counts.len(), work)?;
    row.extend_from_slice(counts);
    output.push(row);
    Ok(())
}
fn vector<T>(length: usize, work: &mut Work) -> Result<Vec<T>> {
    work.consume(length)?;
    work.allocate(
        length
            .checked_mul(size_of::<T>())
            .ok_or_else(|| invalid("allocation size overflow"))?,
    )?;
    let mut value = Vec::new();
    value
        .try_reserve_exact(length)
        .map_err(|_| invalid("allocation failed"))?;
    Ok(value)
}
fn filled<T: Clone>(length: usize, value: T, work: &mut Work) -> Result<Vec<T>> {
    let mut output = vector(length, work)?;
    output.resize(length, value);
    Ok(output)
}
fn count(value: u64) -> Result<u32> {
    u32::try_from(value).map_err(|_| invalid("decomposition multiplicity overflow"))
}
fn increment_count(value: u32) -> Result<u32> {
    value
        .checked_add(1)
        .ok_or_else(|| invalid("decomposition count overflow"))
}
fn product(a: u64, b: u64) -> Result<u64> {
    a.checked_mul(b)
        .ok_or_else(|| invalid("integer product overflow"))
}
fn add(a: u64, b: u64) -> Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| invalid("residue-table mass overflow"))
}
fn gcd(mut a: u64, mut b: u64, work: &mut Work) -> Result<u64> {
    while b != 0 {
        work.consume(1)?;
        (a, b) = (b, a % b);
    }
    Ok(a)
}
fn integer(value: f64) -> Result<u64> {
    if !value.is_finite() || !(0.0..U64_END).contains(&value) {
        Err(invalid("real interval endpoint outside u64 range"))
    } else {
        Ok(value as u64)
    }
}
fn invalid(message: &'static str) -> Error {
    Error::InvalidValue(format!("IMS decomposition: {message}"))
}
struct Work {
    remaining: usize,
    bytes: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_IMS_DECOMPOSER_WORK,
            bytes: MAX_IMS_DECOMPOSER_BYTES,
        }
    }
}
impl Work {
    fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| invalid("work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, amount: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(amount)
            .ok_or_else(|| invalid("logical allocation limit exceeded"))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn shared_budget_failure_does_not_publish_partial_results() {
        let solver = IMSIntegerMassDecomposer::from_integer_weights(&[1, 1]).unwrap();
        let mut work = Work {
            remaining: 40,
            bytes: MAX_IMS_DECOMPOSER_BYTES,
        };
        assert!(solver.all_with_work(100, &mut work).is_err());
        assert!(work.remaining < 5); // Failure after traversing and staging rows.
        assert_eq!(solver.number_of_decompositions(2).unwrap(), 3);
        let mut work = Work {
            remaining: MAX_IMS_DECOMPOSER_WORK,
            bytes: 90,
        };
        assert!(solver.all_with_work(2, &mut work).is_err());
        let real = IMSRealMassDecomposer::new(&IMSWeights::from_masses(&[1.0, 1.0], 1.0).unwrap())
            .unwrap();
        let mut work = Work {
            remaining: 1000,
            bytes: MAX_IMS_DECOMPOSER_BYTES,
        };
        let constraints = BTreeMap::from([(0, (0, 0)), (1, (0, 0))]);
        assert!(
            real.all_with_work(100.0, 1.0, &constraints, &mut work)
                .is_err()
        );
        assert!(work.remaining < 10); // Filtering everything cannot bypass work.
    }
    #[test]
    fn literal_source_cache_counter_can_leave_a_zero_witness() {
        let solver = IMSIntegerMassDecomposer::from_integer_weights(&[10, 6, 15]).unwrap();
        // Source-derived rows: first*last=150, then gcd(10,15)=5. In its
        // second cache loop only counters[1] increments. The update of row[3]
        // from row[8]+15=33 therefore saves untouched counters[3]=0.
        assert_eq!(solver.table[1], [0, 150, 12, 150, 24, 150, 6, 150, 18, 150]);
        assert_eq!(solver.table[2], [0, 21, 12, 33, 24, 15, 6, 27, 18, 39]);
        assert_eq!(solver.witness[3], (2, 0));
        assert!(solver.exists(33).unwrap());
        assert_eq!(solver.decompositions(33).unwrap(), [vec![0, 3, 1]]);
        assert!(solver.decomposition(33).is_err());
        // Nonzero first-coordinate additions do not change the witness residue.
        assert_eq!(solver.decompositions(43).unwrap(), [vec![1, 3, 1]]);
        assert!(solver.decomposition(43).is_err());
    }
    #[test]
    fn construction_and_witness_errors_are_checked() {
        let mut work = Work {
            remaining: 4,
            bytes: MAX_IMS_DECOMPOSER_BYTES,
        };
        assert!(IMSIntegerMassDecomposer::build(&[4, 5, 6], &mut work).is_err());
        let mut work = Work {
            remaining: MAX_IMS_DECOMPOSER_WORK,
            bytes: 8,
        };
        assert!(IMSIntegerMassDecomposer::build(&[4, 5, 6], &mut work).is_err());
        let mut solver = IMSIntegerMassDecomposer::from_integer_weights(&[4, 5, 6]).unwrap();
        solver.witness[2] = (0, 0);
        assert!(solver.decomposition(10).is_err());
        // Finite source break is after witness addition, with no mass check.
        solver.witness[2] = (2, 2);
        assert_eq!(solver.decomposition(10).unwrap(), Some(vec![1, 0, 2]));
        assert!(count(u64::from(u32::MAX) + 1).is_err());
        assert!(integer(U64_END).is_err());
        assert_eq!(
            integer(f64::from_bits(U64_END.to_bits() - 1)).unwrap(),
            u64::MAX - 2047
        );
    }
}
