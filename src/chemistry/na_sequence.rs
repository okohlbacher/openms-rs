// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Owned nucleic-acid sequences with the pinned OpenMS fragment chemistry.
//! Record formulas, rather than their independently declared masses, determine
//! sequence masses. Charges add natural hydrogen and subtract electron mass;
//! the result is an ion mass, not a mass-to-charge ratio.

use super::{
    ELECTRON_MASS_U, EmpiricalFormula, Ribonucleotide, RibonucleotideDB,
    RibonucleotideTermSpecificity,
};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

pub const MAX_NA_SEQUENCE_RESIDUES: usize = 1_000_000;
pub const MAX_NA_SEQUENCE_TEXT_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_NA_SEQUENCE_WORK: usize = 50_000_000;
pub const MAX_NA_SEQUENCE_BYTES: usize = 256 * 1024 * 1024;

/// All finite source fragment types, including its legacy fallback variants.
/// `Internal`, terminal-only, precursor, neutral-loss, and unassigned variants
/// return the backbone formula without end groups or charge hydrogens.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NAFragmentType {
    #[default]
    Full,
    Internal,
    FivePrime,
    ThreePrime,
    AIon,
    BIon,
    CIon,
    XIon,
    YIon,
    ZIon,
    Precursor,
    BIonMinusH2O,
    YIonMinusH2O,
    BIonMinusNH3,
    YIonMinusNH3,
    NonIdentified,
    Unannotated,
    WIon,
    AminusB,
    DIon,
}

/// Immutable record handles; sequence and end replacement are checked and atomic.
/// Equality uses complete chemical record values, never pointer addresses.
#[derive(Clone, Debug)]
pub struct NASequence {
    residues: Vec<Arc<Ribonucleotide>>,
    five_prime: Option<Arc<Ribonucleotide>>,
    three_prime: Option<Arc<Ribonucleotide>>,
    // Resolved once, not looked up globally when slicing a caller-owned sequence.
    phosphorothioate_end: Option<Arc<Ribonucleotide>>,
}

impl Default for NASequence {
    fn default() -> Self {
        Self::new()
    }
}
impl NASequence {
    pub fn new() -> Self {
        Self::empty_with_registry(RibonucleotideDB::global())
    }

    fn empty_with_registry(registry: &RibonucleotideDB) -> Self {
        Self {
            residues: Vec::new(),
            five_prime: None,
            three_prime: None,
            phosphorothioate_end: registry.get("5'-p*").ok(),
        }
    }

    pub fn from_records(residues: Vec<Arc<Ribonucleotide>>) -> Result<Self> {
        Self::from_records_with_registry(residues, RibonucleotideDB::global())
    }

    /// Captures the registry's optional `5'-p*` record for later sulfur slices.
    /// No other record is re-resolved: the supplied owned handles are retained.
    pub fn from_records_with_registry(
        residues: Vec<Arc<Ribonucleotide>>,
        registry: &RibonucleotideDB,
    ) -> Result<Self> {
        let mut result = Self::empty_with_registry(registry);
        result.set_sequence(residues)?;
        Ok(result)
    }

    pub fn parse(input: &str) -> Result<Self> {
        Self::parse_with_registry(input, RibonucleotideDB::global())
    }

    pub fn parse_with_registry(input: &str, registry: &RibonucleotideDB) -> Result<Self> {
        Self::parse_with_work(input, registry, &mut Work::default(), false)
    }

    /// Parse graph parents under one cumulative lookup, work and payload allowance.
    pub(crate) fn parse_with_budget(
        input: &str,
        registry: &RibonucleotideDB,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<Self> {
        let mut work = Work {
            remaining: *remaining_work,
            bytes: *remaining_bytes,
        };
        let result = Self::parse_with_work(input, registry, &mut work, true);
        *remaining_work = work.remaining;
        *remaining_bytes = work.bytes;
        result
    }

    fn parse_with_work(
        input: &str,
        registry: &RibonucleotideDB,
        work: &mut Work,
        charge_payload: bool,
    ) -> Result<Self> {
        if input.len() > MAX_NA_SEQUENCE_TEXT_BYTES {
            return Err(invalid("RNA input exceeds byte limit"));
        }
        work.consume(input.len())?;
        if charge_payload {
            work.parse_payload(std::mem::size_of::<Self>())?;
        }
        let mut result = Self::empty_with_registry(registry);
        if let Some(context) = result
            .phosphorothioate_end
            .as_ref()
            .filter(|_| charge_payload)
        {
            work.parse_payload(context.payload_bytes()?)?;
        }
        if input.is_empty() {
            return Ok(result);
        }
        let mut start = 0;
        let mut stop = input.len();
        match input.as_bytes()[0] {
            b'p' => {
                result.five_prime = Some(lookup(registry, "5'-p", work)?);
                start = 1;
            }
            b'*' => {
                result.five_prime = Some(lookup(registry, "5'-p*", work)?);
                start = 1;
            }
            _ => {}
        }
        if input.len() > 1 {
            match input.as_bytes()[input.len() - 1] {
                b'p' => {
                    result.three_prime = Some(lookup(registry, "3'-p", work)?);
                    stop -= 1;
                }
                b'c' => {
                    result.three_prime = Some(lookup(registry, "3'-c", work)?);
                    stop -= 1;
                }
                _ => {}
            }
        }
        if charge_payload {
            for end in result.five_prime.iter().chain(result.three_prime.iter()) {
                work.parse_payload(end.payload_bytes()?)?;
            }
        }
        while start < stop {
            let byte = input.as_bytes()[start];
            if byte == b' ' {
                start += 1;
                continue;
            }
            let (record, bracketed) = if byte == b'[' {
                let end = input[start + 1..stop]
                    .find(']')
                    .map(|n| start + 1 + n)
                    .ok_or_else(|| invalid("RNA modification is missing ']'"))?;
                let record = lookup(registry, &input[start + 1..end], work)?;
                start = end + 1;
                (record, true)
            } else {
                if !byte.is_ascii() {
                    return Err(invalid("RNA non-ASCII residue codes require brackets"));
                }
                let record = lookup(registry, &input[start..start + 1], work)?;
                start += 1;
                (record, false)
            };
            if charge_payload {
                work.parse_payload(record.payload_bytes()?)?;
            }
            match (bracketed, record.term_specificity()) {
                (true, RibonucleotideTermSpecificity::FivePrime) => {
                    result.five_prime = Some(record)
                }
                (true, RibonucleotideTermSpecificity::ThreePrime) => {
                    result.three_prime = Some(record)
                }
                _ => {
                    if result.residues.len() == MAX_NA_SEQUENCE_RESIDUES {
                        return Err(invalid("RNA residue limit exceeded"));
                    }
                    if charge_payload {
                        // Cover Vec's initial capacity and geometric spare slots
                        // before reserving; records above are charged cumulatively.
                        let slots = if result.residues.is_empty() { 4 } else { 2 };
                        work.parse_payload(slots * std::mem::size_of::<Arc<Ribonucleotide>>())?;
                    }
                    result
                        .residues
                        .try_reserve(1)
                        .map_err(|_| invalid("RNA residue allocation failed"))?;
                    result.residues.push(record);
                }
            }
        }
        if charge_payload {
            work.consume(result.len().saturating_add(3))?;
        }
        result.validate_parts(
            result.residues.iter().map(Arc::as_ref),
            result.len(),
            result.five_prime.as_deref(),
            result.three_prime.as_deref(),
        )?;
        Ok(result)
    }

    pub fn residues(&self) -> &[Arc<Ribonucleotide>] {
        &self.residues
    }
    pub fn len(&self) -> usize {
        self.residues.len()
    }
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }
    /// Clears residues and both attached ends; retains the owned slicing context.
    pub fn clear(&mut self) {
        self.residues.clear();
        self.five_prime = None;
        self.three_prime = None;
    }
    pub fn get_residue(&self, index: usize) -> Result<&Arc<Ribonucleotide>> {
        self.residues
            .get(index)
            .ok_or_else(|| invalid("RNA residue index out of bounds"))
    }
    pub fn set_sequence(&mut self, residues: Vec<Arc<Ribonucleotide>>) -> Result<()> {
        self.validate_parts(
            residues.iter().map(Arc::as_ref),
            residues.len(),
            self.five_prime.as_deref(),
            self.three_prime.as_deref(),
        )?;
        self.residues = residues;
        Ok(())
    }
    pub fn set_residue(&mut self, index: usize, residue: Arc<Ribonucleotide>) -> Result<()> {
        self.get_residue(index)?;
        self.validate_parts(
            self.residues.iter().enumerate().map(|(i, r)| {
                if i == index {
                    residue.as_ref()
                } else {
                    r.as_ref()
                }
            }),
            self.len(),
            self.five_prime.as_deref(),
            self.three_prime.as_deref(),
        )?;
        self.residues[index] = residue;
        Ok(())
    }
    pub fn five_prime_mod(&self) -> Option<&Arc<Ribonucleotide>> {
        self.five_prime.as_ref()
    }
    pub fn three_prime_mod(&self) -> Option<&Arc<Ribonucleotide>> {
        self.three_prime.as_ref()
    }
    /// Source placement is explicit: the record's specificity is not enforced.
    pub fn set_five_prime_mod(&mut self, modification: Option<Arc<Ribonucleotide>>) -> Result<()> {
        self.validate_parts(
            self.residues.iter().map(Arc::as_ref),
            self.len(),
            modification.as_deref(),
            self.three_prime.as_deref(),
        )?;
        self.five_prime = modification;
        Ok(())
    }
    /// Source placement is explicit: the record's specificity is not enforced.
    pub fn set_three_prime_mod(&mut self, modification: Option<Arc<Ribonucleotide>>) -> Result<()> {
        self.validate_parts(
            self.residues.iter().map(Arc::as_ref),
            self.len(),
            self.five_prime.as_deref(),
            modification.as_deref(),
        )?;
        self.three_prime = modification;
        Ok(())
    }

    /// Replaces the owned record used when slicing after a `*` residue.
    /// `None` makes such a slice a checked error; ordinary slices still work.
    pub fn with_phosphorothioate_end(mut self, end: Option<Arc<Ribonucleotide>>) -> Result<Self> {
        if let Some(record) = &end {
            if record.code() != "5'-p*" {
                return Err(invalid("RNA sulfur slicing end must have code 5'-p*"));
            }
        }
        self.phosphorothioate_end = end;
        self.validate_parts(
            self.residues.iter().map(Arc::as_ref),
            self.len(),
            self.five_prime.as_deref(),
            self.three_prime.as_deref(),
        )?;
        Ok(self)
    }

    /// Source rejects a full-length prefix (also zero on an empty sequence).
    pub fn prefix(&self, length: usize) -> Result<Self> {
        if length >= self.len() {
            return Err(invalid(
                "RNA prefix length must be smaller than sequence length",
            ));
        }
        self.slice(0, length, self.five_prime.clone(), None)
    }
    /// Source rejects a full-length suffix; a zero-length suffix can carry ends.
    pub fn suffix(&self, length: usize) -> Result<Self> {
        if length >= self.len() {
            return Err(invalid(
                "RNA suffix length must be smaller than sequence length",
            ));
        }
        let start = self.len() - length;
        self.slice(
            start,
            length,
            self.boundary_end(start)?,
            self.three_prime.clone(),
        )
    }
    /// Oversized length is clamped; `start == len` is an error, even for zero.
    pub fn subsequence(&self, start: usize, length: Option<usize>) -> Result<Self> {
        if start >= self.len() {
            return Err(invalid("RNA subsequence start out of bounds"));
        }
        let length = length.unwrap_or(usize::MAX).min(self.len() - start);
        let five = if start == 0 {
            self.five_prime.clone()
        } else {
            self.boundary_end(start)?
        };
        let three = if start + length == self.len() {
            self.three_prime.clone()
        } else {
            None
        };
        self.slice(start, length, five, three)
    }
    fn boundary_end(&self, start: usize) -> Result<Option<Arc<Ribonucleotide>>> {
        if self.residues[start - 1].code().ends_with('*') {
            self.phosphorothioate_end
                .clone()
                .map(Some)
                .ok_or_else(|| invalid("RNA sulfur slice requires a resolved 5'-p* record"))
        } else {
            Ok(None)
        }
    }
    fn slice(
        &self,
        start: usize,
        length: usize,
        five_prime: Option<Arc<Ribonucleotide>>,
        three_prime: Option<Arc<Ribonucleotide>>,
    ) -> Result<Self> {
        self.validate_parts(
            self.residues[start..start + length].iter().map(Arc::as_ref),
            length,
            five_prime.as_deref(),
            three_prime.as_deref(),
        )?;
        let mut residues = Vec::new();
        residues
            .try_reserve_exact(length)
            .map_err(|_| invalid("RNA slice allocation failed"))?;
        residues.extend_from_slice(&self.residues[start..start + length]);
        Ok(Self {
            residues,
            five_prime,
            three_prime,
            phosphorothioate_end: self.phosphorothioate_end.clone(),
        })
    }

    /// Source formula algebra, including signed atom counts and legacy fallback.
    pub fn formula(&self, fragment: NAFragmentType, charge: i32) -> Result<EmpiricalFormula> {
        self.formula_with_work(fragment, charge, &mut Work::default())
    }
    fn formula_with_work(
        &self,
        fragment: NAFragmentType,
        charge: i32,
        work: &mut Work,
    ) -> Result<EmpiricalFormula> {
        if self.is_empty() {
            return Ok(EmpiricalFormula::default());
        }
        // Preflight the small fixed set of correction formulas and their copies.
        work.consume(128)?;
        work.allocate(64)?;
        let h = EmpiricalFormula::parse("H")?;
        let phosphate = EmpiricalFormula::parse("HPO3")?;
        let thiophosphate = EmpiricalFormula::parse("HPO2S")?;
        let water = EmpiricalFormula::parse("H2O")?;
        let linkage = phosphate.checked_sub(&water)?;
        let sulfur_linkage = thiophosphate.checked_sub(&water)?;
        let mut result = EmpiricalFormula::default();
        for (index, record) in self.residues.iter().enumerate() {
            work.consume(1)?;
            work.combine(&mut result, record.formula(), false)?;
            if index + 1 < self.len() {
                work.combine(
                    &mut result,
                    if record.code().ends_with('*') {
                        &sulfur_linkage
                    } else {
                        &linkage
                    },
                    false,
                )?;
            }
        }
        // The source computes both local end formulas even for fallback types.
        let five = self.local_end(self.five_prime.as_deref(), &h, work)?;
        let three = self.local_end(self.three_prime.as_deref(), &h, work)?;
        use NAFragmentType::*;
        if matches!(
            fragment,
            Internal
                | FivePrime
                | ThreePrime
                | Precursor
                | BIonMinusH2O
                | YIonMinusH2O
                | BIonMinusNH3
                | YIonMinusNH3
                | NonIdentified
                | Unannotated
        ) {
            return Ok(result);
        }
        work.combine(&mut result, &h.checked_scale(charge)?, false)?;
        match fragment {
            Full => {
                work.combine(&mut result, &five, false)?;
                work.combine(&mut result, &three, false)?;
            }
            AminusB => {
                work.combine(&mut result, &five, false)?;
                work.combine(&mut result, &water.checked_scale(-2)?, false)?;
                let last = self.residues.last().expect("nonempty sequence");
                work.combine(&mut result, last.formula(), true)?;
                work.combine(&mut result, last.baseloss_formula(), false)?;
            }
            AIon | BIon | CIon | DIon => {
                work.combine(&mut result, &five, false)?;
                match fragment {
                    AIon => work.combine(&mut result, &water.checked_scale(-1)?, false)?,
                    CIon => {
                        work.combine(&mut result, &EmpiricalFormula::parse("H-1PO2")?, false)?
                    }
                    DIon => work.combine(&mut result, &phosphate, false)?,
                    _ => {}
                }
                if matches!(fragment, CIon | DIon)
                    && self
                        .residues
                        .last()
                        .expect("nonempty sequence")
                        .code()
                        .ends_with('*')
                {
                    work.combine(&mut result, &EmpiricalFormula::parse("SO-1")?, false)?;
                }
            }
            WIon | XIon | YIon | ZIon => {
                work.combine(&mut result, &three, false)?;
                match fragment {
                    WIon => work.combine(&mut result, &phosphate, false)?,
                    XIon => {
                        work.combine(&mut result, &EmpiricalFormula::parse("H-1PO2")?, false)?
                    }
                    ZIon => work.combine(&mut result, &water.checked_scale(-1)?, false)?,
                    _ => {}
                }
                if matches!(fragment, WIon | XIon) && five == thiophosphate {
                    work.combine(&mut result, &EmpiricalFormula::parse("SO-1")?, false)?;
                }
            }
            _ => unreachable!("fallback types returned above"),
        }
        Ok(result)
    }
    fn local_end(
        &self,
        record: Option<&Ribonucleotide>,
        h: &EmpiricalFormula,
        work: &mut Work,
    ) -> Result<EmpiricalFormula> {
        match record {
            Some(record) => {
                work.consume(record.formula().atoms.len())?;
                work.allocate(record.formula().atoms.len())?;
                let mut formula = record.formula().clone();
                work.combine(&mut formula, h, true)?;
                Ok(formula)
            }
            None => Ok(EmpiricalFormula::default()),
        }
    }
    pub fn mono_mass(&self, fragment: NAFragmentType, charge: i32) -> Result<f64> {
        checked_mass(
            self.formula(fragment, charge)?.mono_mass() - f64::from(charge) * ELECTRON_MASS_U,
        )
    }
    /// Return a fragment formula under a caller's cumulative work/allocation allowance.
    pub(crate) fn formula_with_budget(
        &self,
        fragment: NAFragmentType,
        charge: i32,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<EmpiricalFormula> {
        let mut work = Work {
            remaining: *remaining_work,
            bytes: *remaining_bytes,
        };
        let result = self.formula_with_work(fragment, charge, &mut work);
        *remaining_work = work.remaining;
        *remaining_bytes = work.bytes;
        result
    }
    /// Internal callers share their formula work and allocation allowance.
    /// Counters retain work consumed before an error as well as on success.
    pub(crate) fn mono_mass_with_budget(
        &self,
        fragment: NAFragmentType,
        charge: i32,
        remaining_work: &mut usize,
        remaining_bytes: &mut usize,
    ) -> Result<f64> {
        let mut work = Work {
            remaining: *remaining_work,
            bytes: *remaining_bytes,
        };
        let result = self
            .formula_with_work(fragment, charge, &mut work)
            .and_then(|formula| {
                checked_mass(formula.mono_mass() - f64::from(charge) * ELECTRON_MASS_U)
            });
        *remaining_work = work.remaining;
        *remaining_bytes = work.bytes;
        result
    }
    pub fn average_mass(&self, fragment: NAFragmentType, charge: i32) -> Result<f64> {
        checked_mass(
            self.formula(fragment, charge)?.average_mass() - f64::from(charge) * ELECTRON_MASS_U,
        )
    }

    /// Returns source display text only if this registry reconstructs every value
    /// and the sulfur slicing context that can affect this sequence's fragments.
    pub fn checked_string_with_registry(&self, registry: &RibonucleotideDB) -> Result<String> {
        let text = self.to_string();
        let parsed = Self::parse_with_registry(&text, registry)?;
        if &parsed != self {
            return Err(Error::Unsupported("RNA text does not preserve records, placement, or sulfur slicing chemistry in this registry".into()));
        }
        Ok(text)
    }

    /// Logical generation payload, including shared records and slicing context.
    /// Callers charge the traversal before calling and their own output copies.
    pub(crate) fn generation_payload_bytes(&self) -> Result<usize> {
        let initial = self
            .len()
            .checked_mul(std::mem::size_of::<Arc<Ribonucleotide>>())
            .and_then(|bytes| bytes.checked_add(std::mem::size_of::<Self>()))
            .ok_or_else(|| invalid("RNA generation payload overflows"))?;
        self.residues
            .iter()
            .chain(self.five_prime.iter())
            .chain(self.three_prime.iter())
            .chain(self.phosphorothioate_end.iter())
            .try_fold(initial, |bytes, record| {
                bytes
                    .checked_add(record.payload_bytes()?)
                    .ok_or_else(|| invalid("RNA generation payload overflows"))
            })
    }

    /// Validate one final assignment while retaining the caller's slicing record.
    /// The caller preflights vector allocation and charges this validation scan.
    pub(crate) fn with_generation_records(
        &self,
        residues: Vec<Arc<Ribonucleotide>>,
        five_prime: Option<Arc<Ribonucleotide>>,
        three_prime: Option<Arc<Ribonucleotide>>,
    ) -> Result<Self> {
        self.validate_parts(
            residues.iter().map(Arc::as_ref),
            residues.len(),
            five_prime.as_deref(),
            three_prime.as_deref(),
        )?;
        Ok(Self {
            residues,
            five_prime,
            three_prime,
            phosphorothioate_end: self.phosphorothioate_end.clone(),
        })
    }

    fn validate_parts<'a>(
        &self,
        records: impl Iterator<Item = &'a Ribonucleotide>,
        count: usize,
        five: Option<&'a Ribonucleotide>,
        three: Option<&'a Ribonucleotide>,
    ) -> Result<()> {
        if count > MAX_NA_SEQUENCE_RESIDUES {
            return Err(invalid("RNA residue limit exceeded"));
        }
        let mut text = 0usize;
        let mut payload = count
            .checked_mul(std::mem::size_of::<Arc<Ribonucleotide>>())
            .ok_or_else(|| invalid("RNA allocation overflow"))?;
        for record in records.chain(five).chain(three) {
            text = text
                .checked_add(record.code().len() + 2)
                .ok_or_else(|| invalid("RNA text length overflow"))?;
            payload = payload
                .checked_add(record.payload_bytes()?)
                .ok_or_else(|| invalid("RNA payload length overflow"))?;
            if text > MAX_NA_SEQUENCE_TEXT_BYTES || payload > MAX_NA_SEQUENCE_BYTES {
                return Err(invalid("RNA sequence text or payload limit exceeded"));
            }
        }
        if let Some(end) = &self.phosphorothioate_end {
            if payload.saturating_add(end.payload_bytes()?) > MAX_NA_SEQUENCE_BYTES {
                return Err(invalid("RNA sequence payload limit exceeded"));
            }
        }
        Ok(())
    }
    fn relevant_sulfur_end(&self) -> Option<&Arc<Ribonucleotide>> {
        self.residues
            .iter()
            .any(|r| r.code().ends_with('*'))
            .then_some(self.phosphorothioate_end.as_ref())
            .flatten()
    }
}

impl PartialEq for NASequence {
    fn eq(&self, other: &Self) -> bool {
        self.residues == other.residues
            && self.five_prime == other.five_prime
            && self.three_prime == other.three_prime
            && self.relevant_sulfur_end() == other.relevant_sulfur_end()
    }
}
impl Eq for NASequence {}
impl PartialOrd for NASequence {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for NASequence {
    fn cmp(&self, other: &Self) -> Ordering {
        self.five_prime
            .cmp(&other.five_prime)
            .then_with(|| self.len().cmp(&other.len()))
            .then_with(|| {
                self.residues
                    .iter()
                    .map(|r| (r.code(), r.as_ref()))
                    .cmp(other.residues.iter().map(|r| (r.code(), r.as_ref())))
            })
            .then_with(|| self.three_prime.cmp(&other.three_prime))
            .then_with(|| self.relevant_sulfur_end().cmp(&other.relevant_sulfur_end()))
    }
}
impl Hash for NASequence {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.residues.hash(state);
        self.five_prime.hash(state);
        self.three_prime.hash(state);
        self.relevant_sulfur_end().hash(state);
    }
}
impl std::str::FromStr for NASequence {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self> {
        Self::parse(s)
    }
}
impl fmt::Display for NASequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(record) = &self.five_prime {
            match record.code() {
                "5'-p" => f.write_str("p")?,
                "5'-p*" => f.write_str("*")?,
                code => write!(f, "[{code}]")?,
            }
        }
        for record in &self.residues {
            if record.code().len() == 1 {
                f.write_str(record.code())?;
            } else {
                write!(f, "[{}]", record.code())?;
            }
        }
        if let Some(record) = &self.three_prime {
            match record.code() {
                "3'-p" => f.write_str("p")?,
                "3'-c" => f.write_str("c")?,
                code => write!(f, "[{code}]")?,
            }
        }
        Ok(())
    }
}

struct Work {
    remaining: usize,
    bytes: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_NA_SEQUENCE_WORK,
            bytes: MAX_NA_SEQUENCE_BYTES,
        }
    }
}
impl Work {
    fn parse_payload(&mut self, bytes: usize) -> Result<()> {
        self.consume(bytes)?;
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("RNA parent parsing allocation limit exceeded"))?;
        Ok(())
    }
    fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(amount)
            .ok_or_else(|| invalid("RNA operation work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, formula_entries: usize) -> Result<()> {
        // Conservative cumulative BTreeMap node/scratch payload, not atom counts.
        // A sparse root still reserves a whole node, including for one entry.
        self.bytes = self
            .bytes
            .checked_sub(formula_entries.saturating_mul(128).saturating_add(512))
            .ok_or_else(|| invalid("RNA formula allocation limit exceeded"))?;
        Ok(())
    }
    fn combine(
        &mut self,
        left: &mut EmpiricalFormula,
        right: &EmpiricalFormula,
        subtract: bool,
    ) -> Result<()> {
        let n = left.atoms.len().saturating_add(right.atoms.len());
        self.consume(
            n.saturating_mul((usize::BITS - n.max(1).leading_zeros()) as usize + 1)
                .saturating_mul(12),
        )?;
        self.allocate(n)?;
        *left = if subtract {
            left.checked_sub(right)?
        } else {
            left.checked_add(right)?
        };
        Ok(())
    }
}
fn lookup(registry: &RibonucleotideDB, code: &str, work: &mut Work) -> Result<Arc<Ribonucleotide>> {
    let depth = (usize::BITS - registry.entries().len().max(1).leading_zeros()) as usize + 1;
    // A BTreeMap node can compare several keys; 12 comparisons per
    // binary-tree height is a conservative bound for the standard map.
    work.consume((code.len() + 1).saturating_mul(depth).saturating_mul(12))?;
    registry.get(code)
}
fn checked_mass(mass: f64) -> Result<f64> {
    if mass.is_finite() {
        Ok(mass)
    } else {
        Err(invalid("RNA mass is nonfinite"))
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod shared_budget_tests {
    use super::*;

    #[test]
    fn parent_parsing_preflights_allocation_and_shares_counters() {
        let registry = RibonucleotideDB::global();
        let text = "pA[C*]G[3'-p]";
        let expected = NASequence::parse(text).unwrap();
        let (mut work, mut bytes) = (MAX_NA_SEQUENCE_WORK, MAX_NA_SEQUENCE_BYTES);
        let parsed = NASequence::parse_with_budget(text, registry, &mut work, &mut bytes).unwrap();
        assert_eq!(parsed, expected);
        let used_work = MAX_NA_SEQUENCE_WORK - work;
        let used_bytes = MAX_NA_SEQUENCE_BYTES - bytes;
        assert!(used_work > 0 && used_bytes >= parsed.generation_payload_bytes().unwrap());
        for constrain_work in [true, false] {
            let (mut work, mut bytes) = if constrain_work {
                (used_work * 2 - 1, MAX_NA_SEQUENCE_BYTES)
            } else {
                (MAX_NA_SEQUENCE_WORK, used_bytes * 2 - 1)
            };
            assert_eq!(
                NASequence::parse_with_budget(text, registry, &mut work, &mut bytes).unwrap(),
                expected
            );
            assert!(NASequence::parse_with_budget(text, registry, &mut work, &mut bytes).is_err());
        }
        // No sequence payload may be created before a zero byte allowance fails.
        // Exactly the input scan and fixed object preflight have consumed work;
        // no registry lookup, residue allocation or validation was reached.
        let text = "A".repeat(100_000);
        let (mut work, mut bytes) = (MAX_NA_SEQUENCE_WORK, 0);
        assert!(NASequence::parse_with_budget(&text, registry, &mut work, &mut bytes).is_err());
        assert_eq!(
            work,
            MAX_NA_SEQUENCE_WORK - text.len() - std::mem::size_of::<NASequence>()
        );
        assert_eq!(bytes, 0);
    }

    #[test]
    fn mass_calls_share_counters_without_changing_source_arithmetic() {
        let sequence: NASequence = "pA[C*]G[3'-p]".parse().unwrap();
        for fragment in [
            NAFragmentType::Full,
            NAFragmentType::AminusB,
            NAFragmentType::Internal,
        ] {
            for charge in [-3, 0, 2] {
                let (mut work, mut bytes) = (MAX_NA_SEQUENCE_WORK, MAX_NA_SEQUENCE_BYTES);
                let expected = sequence.mono_mass(fragment, charge).unwrap();
                let first = sequence
                    .mono_mass_with_budget(fragment, charge, &mut work, &mut bytes)
                    .unwrap();
                let used_work = MAX_NA_SEQUENCE_WORK - work;
                let used_bytes = MAX_NA_SEQUENCE_BYTES - bytes;
                assert!(used_work > 0 && used_bytes > 0);
                let second = sequence
                    .mono_mass_with_budget(fragment, charge, &mut work, &mut bytes)
                    .unwrap();
                assert_eq!(first.to_bits(), expected.to_bits());
                assert_eq!(second.to_bits(), expected.to_bits());
                assert_eq!(MAX_NA_SEQUENCE_WORK - work, used_work * 2);
                assert_eq!(MAX_NA_SEQUENCE_BYTES - bytes, used_bytes * 2);
            }
        }
    }

    #[test]
    fn errors_preserve_consumed_work_and_allocation_counters() {
        let sequence: NASequence = "A".parse().unwrap();
        let (mut work, mut bytes) = (128, MAX_NA_SEQUENCE_BYTES);
        assert!(
            sequence
                .mono_mass_with_budget(NAFragmentType::Full, 0, &mut work, &mut bytes)
                .is_err()
        );
        assert_eq!(work, 0);
        assert!(bytes < MAX_NA_SEQUENCE_BYTES);
        let (mut work, mut bytes) = (MAX_NA_SEQUENCE_WORK, 0);
        assert!(
            sequence
                .mono_mass_with_budget(NAFragmentType::Full, 0, &mut work, &mut bytes)
                .is_err()
        );
        assert_eq!(work, MAX_NA_SEQUENCE_WORK - 128);
        assert_eq!(bytes, 0);
    }

    #[test]
    fn arithmetic_failure_returns_the_shared_counters_and_leaves_input_unchanged() {
        let record = super::super::RibonucleotideRecord {
            code: "huge".into(),
            formula: "H2147483647".parse().unwrap(),
            ..Default::default()
        };
        let sequence =
            NASequence::from_records(vec![Arc::new(Ribonucleotide::from_record(record).unwrap())])
                .unwrap();
        let original = sequence.clone();
        let (mut work, mut bytes) = (MAX_NA_SEQUENCE_WORK, MAX_NA_SEQUENCE_BYTES);
        assert!(
            sequence
                .mono_mass_with_budget(NAFragmentType::Full, 1, &mut work, &mut bytes)
                .is_err()
        );
        assert!(work < MAX_NA_SEQUENCE_WORK && bytes < MAX_NA_SEQUENCE_BYTES);
        assert_eq!(sequence, original);
    }

    #[test]
    fn empty_source_formula_needs_no_shared_budget() {
        let sequence = NASequence::new();
        let (mut work, mut bytes) = (0, 0);
        let mass = sequence
            .mono_mass_with_budget(NAFragmentType::Full, 2, &mut work, &mut bytes)
            .unwrap();
        assert_eq!(
            mass.to_bits(),
            sequence
                .mono_mass(NAFragmentType::Full, 2)
                .unwrap()
                .to_bits()
        );
        assert_eq!((work, bytes), (0, 0));
    }
}
