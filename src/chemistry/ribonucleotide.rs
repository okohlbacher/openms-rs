// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Complete RNA nucleoside records, following Ribonucleotide.cpp at 7c029e8.
//! Declared masses are independent fields; formula changes never replace them.

use super::EmpiricalFormula;
use crate::concept::hash_utils::{
    fnv1a_hash_bytes, fnv1a_hash_string, hash_char, hash_combine, hash_float, hash_int,
};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};

pub(crate) const MAX_RIBONUCLEOTIDE_CODE_BYTES: usize = 4096;
pub(crate) const MAX_RIBONUCLEOTIDE_TEXT_BYTES: usize = 65_536;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RibonucleotideTermSpecificity {
    #[default]
    Anywhere,
    FivePrime,
    ThreePrime,
}

/// Editable input, frozen and validated by [`Ribonucleotide::from_record`].
#[derive(Clone, Debug, PartialEq)]
pub struct RibonucleotideRecord {
    pub name: String,
    pub code: String,
    pub new_code: String,
    pub html_code: String,
    pub formula: EmpiricalFormula,
    pub origin: char,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub term_specificity: RibonucleotideTermSpecificity,
    pub baseloss_formula: EmpiricalFormula,
}

impl Default for RibonucleotideRecord {
    fn default() -> Self {
        Self {
            name: "unknown ribonucleotide".into(),
            code: ".".into(),
            new_code: String::new(),
            html_code: ".".into(),
            formula: EmpiricalFormula::default(),
            origin: '.',
            mono_mass: 0.0,
            average_mass: 0.0,
            term_specificity: RibonucleotideTermSpecificity::Anywhere,
            baseloss_formula: "C5H10O5".parse().expect("constant base-loss formula"),
        }
    }
}

/// Immutable, fully comparable RNA record. Finite signed masses are retained.
/// A code must be nonempty and at most 4096 bytes; each other text field is
/// limited to 65,536 bytes. These bounds also protect registry/sequence copying.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Ribonucleotide {
    record: RibonucleotideRecord,
}

impl Ribonucleotide {
    pub fn from_record(record: RibonucleotideRecord) -> Result<Self> {
        if record.code.is_empty() || record.code.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES {
            return Err(invalid("ribonucleotide code must contain 1..=4096 bytes"));
        }
        if [&record.name, &record.new_code, &record.html_code]
            .iter()
            .any(|text| text.len() > MAX_RIBONUCLEOTIDE_TEXT_BYTES)
        {
            return Err(invalid("ribonucleotide text field exceeds 65536 bytes"));
        }
        if !record.mono_mass.is_finite() || !record.average_mass.is_finite() {
            return Err(invalid("ribonucleotide declared masses must be finite"));
        }
        Ok(Self { record })
    }

    pub fn name(&self) -> &str {
        &self.record.name
    }
    pub fn code(&self) -> &str {
        &self.record.code
    }
    pub fn new_code(&self) -> &str {
        &self.record.new_code
    }
    pub fn html_code(&self) -> &str {
        &self.record.html_code
    }
    pub fn formula(&self) -> &EmpiricalFormula {
        &self.record.formula
    }
    pub fn origin(&self) -> char {
        self.record.origin
    }
    pub fn mono_mass(&self) -> f64 {
        self.record.mono_mass
    }
    pub fn average_mass(&self) -> f64 {
        self.record.average_mass
    }
    pub fn term_specificity(&self) -> RibonucleotideTermSpecificity {
        self.record.term_specificity
    }
    pub fn baseloss_formula(&self) -> &EmpiricalFormula {
        &self.record.baseloss_formula
    }

    /// Source byte-length test; a Unicode multibyte code is considered modified.
    pub fn is_modified(&self) -> bool {
        self.code().len() != 1 || !self.code().starts_with(self.origin())
    }

    /// Source record predicate recognizes `?`, but deliberately not `?*`.
    pub fn is_ambiguous(&self) -> bool {
        self.code().ends_with('?')
    }

    /// Conservative logical payload, including strings, formula nodes and Arc.
    pub(crate) fn payload_bytes(&self) -> Result<usize> {
        [
            self.name().len(),
            self.code().len(),
            self.new_code().len(),
            self.html_code().len(),
        ]
        .into_iter()
        .try_fold(std::mem::size_of::<Self>() + 16, |sum, bytes| {
            sum.checked_add(bytes)
        })
        .and_then(|bytes| {
            self.formula()
                .atoms
                .len()
                .checked_add(self.baseloss_formula().atoms.len())
                .and_then(|n| n.checked_mul(64))
                .and_then(|n| bytes.checked_add(n))
        })
        .ok_or_else(|| invalid("ribonucleotide payload size overflows"))
    }
}

impl Ribonucleotide {
    /// The text source `operator<<(std::ostream&, const Ribonucleotide&)`
    /// writes: `Ribonucleotide '<code>' (<name>, <formula>)`, the formula in
    /// the source's stream form ([`EmpiricalFormula::source_stream_text`], which
    /// omits counts of one).
    ///
    /// [`std::fmt::Display`] writes the same frame with the formula's native
    /// round-trip form instead, so `C10H12N5O6P` there reads `C10H12N5O6P1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::chemistry::ribonucleotide::{Ribonucleotide, RibonucleotideRecord};
    ///
    /// let record = RibonucleotideRecord {
    ///     name: "adenosine".into(),
    ///     code: "A".into(),
    ///     formula: "C10H13N5O4".parse().unwrap(),
    ///     ..RibonucleotideRecord::default()
    /// };
    /// let adenosine = Ribonucleotide::from_record(record).unwrap();
    /// assert_eq!(
    ///     adenosine.source_stream_text(),
    ///     "Ribonucleotide 'A' (adenosine, C10H13N5O4)"
    /// );
    /// ```
    pub fn source_stream_text(&self) -> String {
        format!(
            "Ribonucleotide '{}' ({}, {})",
            self.code(),
            self.name(),
            self.formula().source_stream_text()
        )
    }

    /// The value source `std::hash<Ribonucleotide>` computes.
    ///
    /// Starting from a zero seed, [`hash_combine`] folds in, in order, the
    /// [`fnv1a_hash_string`] of name, code, new code and HTML code, the
    /// [`EmpiricalFormula::source_hash64`] of the formula and of the base-loss
    /// formula, [`hash_char`] of the origin, [`hash_float`] of the mono and
    /// average masses, and [`hash_int`] of the terminal specificity's source
    /// enum value (`ANYWHERE` 0, `FIVE_PRIME` 1, `THREE_PRIME` 2) as a 4-byte
    /// `int`.
    ///
    /// The source origin is one `char` byte. An origin above U+00FF has no
    /// single-byte form; its UTF-8 bytes are hashed with [`fnv1a_hash_bytes`]
    /// instead, which a source record cannot produce.
    pub fn source_hash64(&self) -> u64 {
        let mut seed = 0u64;
        for text in [self.name(), self.code(), self.new_code(), self.html_code()] {
            hash_combine(&mut seed, fnv1a_hash_string(text));
        }
        hash_combine(&mut seed, self.formula().source_hash64());
        hash_combine(&mut seed, self.baseloss_formula().source_hash64());
        let origin = match u8::try_from(u32::from(self.origin())) {
            Ok(byte) => hash_char(byte),
            Err(_) => fnv1a_hash_bytes(self.origin().encode_utf8(&mut [0; 4]).as_bytes()),
        };
        hash_combine(&mut seed, origin);
        hash_combine(&mut seed, hash_float(self.mono_mass()));
        hash_combine(&mut seed, hash_float(self.average_mass()));
        let term: i32 = match self.term_specificity() {
            RibonucleotideTermSpecificity::Anywhere => 0,
            RibonucleotideTermSpecificity::FivePrime => 1,
            RibonucleotideTermSpecificity::ThreePrime => 2,
        };
        hash_combine(&mut seed, hash_int(term));
        seed
    }
}

impl Eq for Ribonucleotide {}
impl Ord for Ribonucleotide {
    fn cmp(&self, other: &Self) -> Ordering {
        self.name()
            .cmp(other.name())
            .then_with(|| self.code().cmp(other.code()))
            .then_with(|| self.new_code().cmp(other.new_code()))
            .then_with(|| self.html_code().cmp(other.html_code()))
            .then_with(|| formula_order(self.formula(), other.formula()))
            .then_with(|| self.origin().cmp(&other.origin()))
            .then_with(|| mass_order(self.mono_mass(), other.mono_mass()))
            .then_with(|| mass_order(self.average_mass(), other.average_mass()))
            .then_with(|| self.term_specificity().cmp(&other.term_specificity()))
            .then_with(|| formula_order(self.baseloss_formula(), other.baseloss_formula()))
    }
}
impl PartialOrd for Ribonucleotide {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Hash for Ribonucleotide {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.name().hash(state);
        self.code().hash(state);
        self.new_code().hash(state);
        self.html_code().hash(state);
        hash_formula(self.formula(), state);
        self.origin().hash(state);
        mass_bits(self.mono_mass()).hash(state);
        mass_bits(self.average_mass()).hash(state);
        self.term_specificity().hash(state);
        hash_formula(self.baseloss_formula(), state);
    }
}
impl fmt::Display for Ribonucleotide {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            output,
            "Ribonucleotide '{}' ({}, {})",
            self.code(),
            self.name(),
            self.formula()
        )
    }
}
fn formula_order(left: &EmpiricalFormula, right: &EmpiricalFormula) -> Ordering {
    left.atoms
        .cmp(&right.atoms)
        .then_with(|| left.charge.cmp(&right.charge))
}
fn hash_formula<H: Hasher>(formula: &EmpiricalFormula, state: &mut H) {
    formula.atoms.len().hash(state);
    for (atom, count) in &formula.atoms {
        atom.symbol.hash(state);
        atom.isotope.hash(state);
        count.hash(state);
    }
    formula.charge.hash(state);
}
fn mass_bits(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}
fn mass_order(left: f64, right: f64) -> Ordering {
    if left == right {
        Ordering::Equal
    } else {
        left.total_cmp(&right)
    }
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
