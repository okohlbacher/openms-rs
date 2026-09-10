// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Modified-RNA digestion using the fourteen pinned enzyme definitions.
//! Per-code predicates preserve source regex search semantics. Arbitrary Boost
//! regex syntax remains outside the native subset.

use super::{NASequence, Ribonucleotide, RibonucleotideDB};
use crate::identification::graph::{
    IdentificationData, IdentifiedOligo, MoleculeType, ParentMatch,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::mem::size_of;
use std::sync::{Arc, OnceLock};

const MAX_ENZYME_RECORDS: usize = 100_000;
const MAX_ENZYME_BYTES: usize = 128 * 1024 * 1024;
const MAX_TEXT_BYTES: usize = 65_536;
const MAX_PATTERNS: usize = 1_024;
pub const MAX_RNASE_WORK: usize = 50_000_000;
pub const MAX_RNASE_PRODUCTS: usize = 100_000;
pub const MAX_RNASE_OUTPUT_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_RNASE_RESIDUES: usize = 1_000_000;

/// Editable source fields; constructing an enzyme freezes their complete values.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DigestionEnzymeRNARecord {
    pub name: String,
    pub regex: String,
    pub regex_description: String,
    pub synonyms: BTreeSet<String>,
    pub cuts_after: String,
    pub cuts_before: String,
    pub five_prime_gain: String,
    pub three_prime_gain: String,
}
impl Default for DigestionEnzymeRNARecord {
    fn default() -> Self {
        Self {
            name: "unknown_enzyme".into(),
            regex: String::new(),
            regex_description: String::new(),
            synonyms: BTreeSet::new(),
            cuts_after: String::new(),
            cuts_before: String::new(),
            five_prime_gain: String::new(),
            three_prime_gain: String::new(),
        }
    }
}

/// RNA fields participate in native value identity, unlike the inherited C++ equality.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DigestionEnzymeRNA {
    record: DigestionEnzymeRNARecord,
}
impl DigestionEnzymeRNA {
    pub fn from_record(record: DigestionEnzymeRNARecord) -> Result<Self> {
        let enzyme = Self { record };
        enzyme.payload_bytes()?;
        Ok(enzyme)
    }
    pub fn record(&self) -> &DigestionEnzymeRNARecord {
        &self.record
    }
    pub fn name(&self) -> &str {
        &self.record.name
    }
    /// Inherited main regex; distinct from modification-aware per-code patterns.
    pub fn regex(&self) -> &str {
        &self.record.regex
    }
    pub fn regex_description(&self) -> &str {
        &self.record.regex_description
    }
    pub fn synonyms(&self) -> &BTreeSet<String> {
        &self.record.synonyms
    }
    pub fn cuts_after(&self) -> &str {
        &self.record.cuts_after
    }
    pub fn cuts_before(&self) -> &str {
        &self.record.cuts_before
    }
    pub fn five_prime_gain(&self) -> &str {
        &self.record.five_prime_gain
    }
    pub fn three_prime_gain(&self) -> &str {
        &self.record.three_prime_gain
    }
    fn payload_bytes(&self) -> Result<usize> {
        if self.synonyms().len() > MAX_PATTERNS {
            return Err(invalid("too many enzyme synonyms"));
        }
        [
            self.name(),
            self.regex(),
            self.regex_description(),
            self.cuts_after(),
            self.cuts_before(),
            self.five_prime_gain(),
            self.three_prime_gain(),
        ]
        .into_iter()
        .chain(self.synonyms().iter().map(String::as_str))
        .try_fold(size_of::<Self>() + self.synonyms().len() * 64, |sum, s| {
            if s.len() > MAX_TEXT_BYTES {
                return Err(invalid("enzyme text exceeds byte limit"));
            }
            sum.checked_add(s.len())
                .ok_or_else(|| invalid("enzyme payload overflow"))
        })
    }
}

/// Owned registry with later-provider name replacement and source alias rules.
/// Iteration uses stable surviving provider order instead of pointer addresses.
#[derive(Clone, Debug, Default)]
pub struct RNaseDB {
    enzymes: Vec<Arc<DigestionEnzymeRNA>>,
    names: BTreeMap<String, Arc<DigestionEnzymeRNA>>,
    regexes: BTreeMap<String, Arc<DigestionEnzymeRNA>>,
}
impl RNaseDB {
    pub fn from_records(records: Vec<DigestionEnzymeRNA>) -> Result<Self> {
        if records.len() > MAX_ENZYME_RECORDS {
            return Err(invalid("too many enzymes"));
        }
        let mut result = Self::default();
        let mut bytes = 0usize;
        let mut budget = Budget(MAX_RNASE_WORK);
        for record in records {
            let payload = record.payload_bytes()?;
            budget.charge(payload)?;
            // Count all input records and duplicate index strings conservatively.
            bytes = add(
                bytes,
                payload
                    .checked_mul(4)
                    .ok_or_else(|| invalid("enzyme payload overflow"))?,
            )?;
            if bytes > MAX_ENZYME_BYTES {
                return Err(invalid("enzyme registry byte limit exceeded"));
            }
            let name = record.name();
            budget.tree_key(name, result.names.len())?;
            if let Some(old) = result.names.get(name).cloned() {
                budget.charge(result.enzymes.len())?;
                result.enzymes.retain(|e| !Arc::ptr_eq(e, &old));
                budget.tree_key(old.name(), result.names.len())?;
                result.names.remove(old.name());
                budget.tree_key(old.name(), result.names.len())?;
                result.names.remove(&old.name().to_ascii_lowercase());
                for alias in old.synonyms() {
                    budget.tree_key(alias, result.names.len())?;
                    result.names.remove(alias);
                }
                if !old.regex().is_empty() {
                    budget.tree_key(old.regex(), result.regexes.len())?;
                    result.regexes.remove(old.regex());
                }
            }
            let record = Arc::new(record);
            budget.tree_key(record.name(), result.names.len())?;
            result.names.insert(record.name().into(), record.clone());
            budget.tree_key(record.name(), result.names.len())?;
            result
                .names
                .insert(record.name().to_ascii_lowercase(), record.clone());
            for alias in record.synonyms() {
                budget.tree_key(alias, result.names.len())?;
                result.names.insert(alias.clone(), record.clone());
            }
            if !record.regex().is_empty() {
                budget.tree_key(record.regex(), result.regexes.len())?;
                result.regexes.insert(record.regex().into(), record.clone());
            }
            result.enzymes.push(record);
        }
        Ok(result)
    }
    pub fn global() -> &'static Self {
        static DB: OnceLock<RNaseDB> = OnceLock::new();
        DB.get_or_init(|| {
            let data: &[(&str, &str, &str, &str, &str, &str)] =
                include!("../../resources/enzymes/rna_enzymes.rs");
            let records = data
                .iter()
                .map(|&(name, description, after, before, five, three)| {
                    DigestionEnzymeRNA::from_record(DigestionEnzymeRNARecord {
                        name: name.into(),
                        regex_description: description.into(),
                        cuts_after: after.into(),
                        cuts_before: before.into(),
                        five_prime_gain: five.into(),
                        three_prime_gain: three.into(),
                        ..DigestionEnzymeRNARecord::default()
                    })
                    .expect("validated pinned RNA enzyme")
                })
                .collect();
            Self::from_records(records).expect("validated pinned RNA enzyme registry")
        })
    }
    pub fn enzymes(&self) -> &[Arc<DigestionEnzymeRNA>] {
        &self.enzymes
    }
    pub fn get_enzyme(&self, name: &str) -> Result<Arc<DigestionEnzymeRNA>> {
        if name.len() > MAX_TEXT_BYTES {
            return Err(invalid("enzyme lookup exceeds byte limit"));
        }
        self.names
            .get(name)
            .cloned()
            .ok_or_else(|| invalid("unknown RNA enzyme name"))
    }
    pub fn has_enzyme(&self, name: &str) -> bool {
        self.names.contains_key(name)
    }
    pub fn names(&self) -> Vec<&str> {
        self.enzymes.iter().map(|e| e.name()).collect()
    }
    pub fn has_regex(&self, regex: &str) -> bool {
        self.regexes.contains_key(regex)
    }
    pub fn enzyme_by_regex(&self, regex: &str) -> Result<Arc<DigestionEnzymeRNA>> {
        self.regexes
            .get(regex)
            .cloned()
            .ok_or_else(|| invalid("unregistered RNA main regex"))
    }
}

#[derive(Clone, Debug)]
enum Pattern {
    Any,
    NoC,
    NotM(&'static [u8]),
    GOrNotM(u8),
    AnyOf(&'static [u8]),
    Rnase4p,
    NotPreceded(&'static [u8], u8),
    Literal(String),
}
impl Pattern {
    fn compile(text: &str) -> Result<Self> {
        if text.len() > MAX_TEXT_BYTES {
            return Err(invalid("RNA pattern exceeds byte limit"));
        }
        Ok(match text {
            "" | ".*" | ".*(?!m)$" => Self::Any,
            "^[^C]+$" => Self::NoC,
            "G(?!m)" => Self::NotM(b"G"),
            "C(?!m)" => Self::NotM(b"C"),
            "[CUY](?!m)" => Self::NotM(b"CUY"),
            "G|A(?!m)" => Self::GOrNotM(b'A'),
            "G|Q(?!m)" => Self::GOrNotM(b'Q'),
            "U|Y" => Self::AnyOf(b"UY"),
            "G|A" => Self::AnyOf(b"GA"),
            "U|P|\\]|D|5" => Self::AnyOf(b"UP]D5"),
            "U|P|m1Y|D|5" => Self::Rnase4p,
            "(?<!m6)A" => Self::NotPreceded(b"m6", b'A'),
            "(?<!m5)C" => Self::NotPreceded(b"m5", b'C'),
            _ if !text.bytes().any(|b| b"\\.^$[]()|*+?{}".contains(&b)) => {
                Self::Literal(text.into())
            }
            _ => {
                return Err(Error::Unsupported(
                    "RNA code pattern is outside the pinned expressions and literal-search subset"
                        .into(),
                ));
            }
        })
    }
    fn matches(&self, code: &str, budget: &mut Budget) -> Result<bool> {
        let bytes = code.as_bytes();
        let width = if let Self::Literal(s) = self {
            s.len().max(1)
        } else {
            4
        };
        budget.charge(
            code.len()
                .checked_mul(width)
                .and_then(|n| n.checked_add(1))
                .ok_or_else(|| invalid("RNA predicate work overflow"))?,
        )?;
        let not_m = |letters: &[u8]| {
            bytes
                .iter()
                .enumerate()
                .any(|(i, b)| letters.contains(b) && bytes.get(i + 1) != Some(&b'm'))
        };
        Ok(match self {
            Self::Any => true,
            Self::NoC => !bytes.is_empty() && !bytes.contains(&b'C'),
            Self::NotM(letters) => not_m(letters),
            Self::GOrNotM(letter) => bytes.contains(&b'G') || not_m(&[*letter]),
            Self::AnyOf(letters) => bytes.iter().any(|b| letters.contains(b)),
            Self::Rnase4p => bytes.iter().any(|b| b"UPD5".contains(b)) || code.contains("m1Y"),
            Self::NotPreceded(prefix, letter) => bytes.iter().enumerate().any(|(i, b)| {
                b == letter && (i < prefix.len() || &bytes[i - prefix.len()..i] != *prefix)
            }),
            Self::Literal(text) => code.contains(text),
        })
    }
}

/// Match one registered per-code expression or a literal substring, without origin substitution.
pub fn matches_code_pattern(pattern: &str, code: &str) -> Result<bool> {
    if code.len() > MAX_TEXT_BYTES {
        return Err(invalid("RNA code exceeds byte limit"));
    }
    Pattern::compile(pattern)?.matches(code, &mut Budget(MAX_RNASE_WORK))
}
fn compile_patterns(text: &str) -> Result<Vec<Pattern>> {
    if text.is_empty() {
        return Ok(Vec::new());
    }
    if text.split(',').count() > MAX_PATTERNS {
        return Err(invalid("too many RNA code patterns"));
    }
    text.split(',').map(Pattern::compile).collect()
}

/// A digestion product with a zero-based half-open source residue interval.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DigestedOligo {
    pub sequence: NASequence,
    pub start: usize,
    pub end: usize,
    pub missed_cleavages: usize,
}

/// Modification-aware RNA digestion; defaults to a configured RNase_T1.
#[derive(Clone, Debug)]
pub struct RNaseDigestion {
    enzyme: Arc<DigestionEnzymeRNA>,
    five: Option<Arc<Ribonucleotide>>,
    three: Option<Arc<Ribonucleotide>>,
    after: Vec<Pattern>,
    before: Vec<Pattern>,
    pub missed_cleavages: usize,
    /// Zero means one residue, matching the source.
    pub min_length: usize,
    /// Zero means the full input length.
    pub max_length: usize,
    pub max_residues: usize,
    pub max_products: usize,
    pub max_work: usize,
    pub max_output_bytes: usize,
}
impl Default for RNaseDigestion {
    fn default() -> Self {
        Self::new("RNase_T1").expect("validated pinned default RNase")
    }
}
impl RNaseDigestion {
    pub fn new(name: &str) -> Result<Self> {
        Self::with_enzyme(
            RNaseDB::global().get_enzyme(name)?,
            RibonucleotideDB::global(),
        )
    }
    pub fn with_enzyme(
        enzyme: Arc<DigestionEnzymeRNA>,
        registry: &RibonucleotideDB,
    ) -> Result<Self> {
        let (five, three, after, before) = configure(&enzyme, registry)?;
        Ok(Self {
            enzyme,
            five,
            three,
            after,
            before,
            missed_cleavages: 0,
            min_length: 0,
            max_length: 0,
            max_residues: MAX_RNASE_RESIDUES,
            max_products: MAX_RNASE_PRODUCTS,
            max_work: MAX_RNASE_WORK,
            max_output_bytes: MAX_RNASE_OUTPUT_BYTES,
        })
    }
    pub fn enzyme(&self) -> &Arc<DigestionEnzymeRNA> {
        &self.enzyme
    }
    pub fn set_enzyme(&mut self, name: &str) -> Result<()> {
        self.set_enzyme_with_registry(
            RNaseDB::global().get_enzyme(name)?,
            RibonucleotideDB::global(),
        )
    }
    pub fn set_enzyme_with_registry(
        &mut self,
        enzyme: Arc<DigestionEnzymeRNA>,
        registry: &RibonucleotideDB,
    ) -> Result<()> {
        let (five, three, after, before) = configure(&enzyme, registry)?;
        self.enzyme = enzyme;
        self.five = five;
        self.three = three;
        self.after = after;
        self.before = before;
        Ok(())
    }
    pub fn digest(&self, sequence: &NASequence) -> Result<Vec<NASequence>> {
        let products = self.digest_with_positions(sequence)?;
        let mut sequences = Vec::new();
        sequences
            .try_reserve_exact(products.len())
            .map_err(|_| invalid("RNA sequence output allocation failed"))?;
        sequences.extend(products.into_iter().map(|p| p.sequence));
        Ok(sequences)
    }
    /// Replace only after all products pass their chemistry and resource checks.
    pub fn digest_into(&self, sequence: &NASequence, output: &mut Vec<NASequence>) -> Result<()> {
        *output = self.digest(sequence)?;
        Ok(())
    }
    pub fn digest_with_positions(&self, sequence: &NASequence) -> Result<Vec<DigestedOligo>> {
        let mut budget = Budget(self.max_work);
        self.digest_with_resources(sequence, &mut budget, &mut 0, &mut 0)
    }

    /// Register all RNA parents' products, retaining existing oligos and processing history.
    /// Equal chemical products merge their inclusive parent positions. The entire
    /// operation is atomic and shares parsing/digestion limits across all parents.
    pub fn digest_identification_data(&self, graph: &mut IdentificationData) -> Result<()> {
        self.digest_identification_data_with_registry(graph, RibonucleotideDB::global())
    }

    /// Registry-aware graph digestion preserves caller-owned modified RNA chemistry.
    pub fn digest_identification_data_with_registry(
        &self,
        graph: &mut IdentificationData,
        registry: &RibonucleotideDB,
    ) -> Result<()> {
        self.validate_limits(0)?;
        let mut budget = Budget(self.max_work);
        budget.charge(graph.parent_count())?;
        let mut parents = Vec::new();
        let mut retained_bytes = 0usize;
        for (id, parent) in graph.parents() {
            if parent.molecule_type == MoleculeType::RNA {
                retained_bytes = add(retained_bytes, size_of_val(&id))?;
                if retained_bytes > self.max_output_bytes {
                    return Err(invalid("RNA parent list byte limit exceeded"));
                }
                parents
                    .try_reserve(1)
                    .map_err(|_| invalid("RNA parent list allocation failed"))?;
                parents.push(id);
            }
        }
        let mut residues = 0usize;
        let mut products = 0usize;
        graph.transaction(|staged| {
            for parent_id in parents {
                // Parsing and sequence generation use the same remaining work;
                // parsed parent payload also consumes the cumulative byte allowance.
                let mut remaining_bytes = self.max_output_bytes - retained_bytes;
                let sequence = NASequence::parse_with_budget(
                    &staged.parent(parent_id)?.sequence,
                    registry,
                    &mut budget.0,
                    &mut remaining_bytes,
                )?;
                retained_bytes = self.max_output_bytes - remaining_bytes;
                residues = add(residues, sequence.len())?;
                self.validate_limits(residues)?;
                for product in self.digest_with_resources(
                    &sequence,
                    &mut budget,
                    &mut retained_bytes,
                    &mut products,
                )? {
                    budget.charge(1)?;
                    let end = product
                        .end
                        .checked_sub(1)
                        .ok_or_else(|| invalid("empty RNA graph digestion product"))?;
                    let mut parent_match = ParentMatch::new(Some(product.start), Some(end));
                    parent_match.left_neighbor = if product.start == 0 {
                        "[".into()
                    } else {
                        neighbor_code(&sequence.residues()[product.start - 1])?
                    };
                    parent_match.right_neighbor = if product.end == sequence.len() {
                        "]".into()
                    } else {
                        neighbor_code(&sequence.residues()[product.end])?
                    };
                    let mut oligo = IdentifiedOligo::new(product.sequence);
                    oligo
                        .parent_matches
                        .insert(parent_id, BTreeSet::from([parent_match]));
                    staged.register_identified_oligo(oligo)?;
                }
            }
            Ok(())
        })
    }
    fn digest_with_resources(
        &self,
        sequence: &NASequence,
        budget: &mut Budget,
        retained_bytes: &mut usize,
        product_count: &mut usize,
    ) -> Result<Vec<DigestedOligo>> {
        self.validate_limits(sequence.len())?;
        budget.charge(sequence.len())?;
        let remaining_products = self
            .max_products
            .checked_sub(*product_count)
            .ok_or_else(|| invalid("RNA digestion product limit exceeded"))?;
        let positions = self.positions(sequence, budget, remaining_products)?;
        if positions.is_empty() {
            return Ok(Vec::new());
        }
        // Prefix accounting measures actual fragment residues, avoiding a full-input
        // payload estimate per tiny product. End/context allowance remains conservative.
        let mut prefix = Vec::new();
        prefix
            .try_reserve_exact(sequence.len() + 1)
            .map_err(|_| invalid("RNA prefix allocation failed"))?;
        prefix.push(0usize);
        for record in sequence.residues() {
            let bytes = add(record.payload_bytes()?, size_of::<Arc<Ribonucleotide>>())?;
            budget.charge(bytes)?;
            prefix.push(add(*prefix.last().expect("prefix origin"), bytes)?);
        }
        budget.charge(sequence.len())?;
        let overhead = sequence
            .generation_payload_bytes()?
            .checked_sub(*prefix.last().expect("prefix total"))
            .ok_or_else(|| invalid("RNA payload accounting mismatch"))?;
        let gains = self
            .five
            .iter()
            .chain(self.three.iter())
            .try_fold(0usize, |n, r| add(n, r.payload_bytes()?))?;
        for &(start, end, _) in &positions {
            budget.charge(end - start + 1)?;
            *retained_bytes = add(
                *retained_bytes,
                add(
                    add(prefix[end] - prefix[start], overhead)?,
                    add(gains, size_of::<DigestedOligo>())?,
                )?,
            )?;
            if *retained_bytes > self.max_output_bytes {
                return Err(invalid("RNA digestion output byte limit exceeded"));
            }
        }
        *product_count = add(*product_count, positions.len())?;
        let mut result = Vec::new();
        result
            .try_reserve_exact(positions.len())
            .map_err(|_| invalid("RNA digestion allocation failed"))?;
        for (start, end, missed_cleavages) in positions {
            // Slicing must occur before gain replacement: a missing sulfur context
            // is a source slicing error even if the enzyme would overwrite that end.
            budget.charge(
                (end - start)
                    .checked_mul(4)
                    .ok_or_else(|| invalid("RNA copy work overflow"))?,
            )?;
            let mut fragment = sequence.subsequence(start, Some(end - start))?;
            if start > 0 {
                fragment.set_five_prime_mod(self.five.clone())?;
            }
            if end < sequence.len() {
                fragment.set_three_prime_mod(self.three.clone())?;
            }
            result.push(DigestedOligo {
                sequence: fragment,
                start,
                end,
                missed_cleavages,
            });
        }
        Ok(result)
    }
    fn validate_limits(&self, length: usize) -> Result<()> {
        if self.max_residues > MAX_RNASE_RESIDUES
            || length > self.max_residues
            || self.max_products > MAX_RNASE_PRODUCTS
            || self.max_output_bytes > MAX_RNASE_OUTPUT_BYTES
            || self.max_work > MAX_RNASE_WORK
        {
            return Err(invalid("RNA digestion resource limits exceeded"));
        }
        Ok(())
    }
    fn positions(
        &self,
        sequence: &NASequence,
        budget: &mut Budget,
        max_products: usize,
    ) -> Result<Vec<(usize, usize, usize)>> {
        let n = sequence.len();
        let minimum = self.min_length.max(1);
        let maximum = if self.max_length == 0 {
            n
        } else {
            self.max_length.min(n)
        };
        let mut result = Vec::new();
        if n == 0 || minimum > maximum {
            return Ok(result);
        }
        let mut emit =
            |start: usize, end: usize, missed: usize, budget: &mut Budget| -> Result<()> {
                budget.charge(1)?;
                if (minimum..=maximum).contains(&(end - start)) {
                    if result.len() >= max_products {
                        return Err(invalid("RNA digestion product limit exceeded"));
                    }
                    if (result.len() + 1).saturating_mul(size_of::<(usize, usize, usize)>())
                        > self.max_output_bytes
                    {
                        return Err(invalid("RNA digestion plan byte limit exceeded"));
                    }
                    result
                        .try_reserve(1)
                        .map_err(|_| invalid("RNA digestion plan allocation failed"))?;
                    result.push((start, end, missed));
                }
                Ok(())
            };
        match self.enzyme.name() {
            "no cleavage" => emit(0, n, 0, budget)?,
            "unspecific cleavage" => {
                for start in 0..=n - minimum {
                    for end in start + minimum..=start + maximum.min(n - start) {
                        emit(start, end, 0, budget)?;
                    }
                }
            }
            _ => {
                let mut cuts = Vec::new();
                cuts.try_reserve_exact(n + 1)
                    .map_err(|_| invalid("RNA cut allocation failed"))?;
                cuts.push(0);
                for i in 1..n {
                    budget.charge(1)?;
                    if i < self.after.len() || n - i < self.before.len() {
                        continue;
                    }
                    let mut matched = true;
                    for (offset, pattern) in self.after.iter().enumerate() {
                        if !pattern.matches(
                            sequence.residues()[i - self.after.len() + offset].code(),
                            budget,
                        )? {
                            matched = false;
                            break;
                        }
                    }
                    if matched {
                        for (offset, pattern) in self.before.iter().enumerate() {
                            if !pattern.matches(sequence.residues()[i + offset].code(), budget)? {
                                matched = false;
                                break;
                            }
                        }
                    }
                    if matched {
                        cuts.push(i);
                    }
                }
                cuts.push(n);
                // Source is start-first, then missed cleavages; protein digestion differs.
                for start in 0..cuts.len() - 1 {
                    let max_missed = self.missed_cleavages.min(cuts.len() - start - 2);
                    for missed in 0..=max_missed {
                        // Later endpoints only grow; rejected long products cannot
                        // reenter the requested length range.
                        if cuts[start + missed + 1] - cuts[start] > maximum {
                            break;
                        }
                        emit(cuts[start], cuts[start + missed + 1], missed, budget)?;
                    }
                }
            }
        }
        Ok(result)
    }
}

fn neighbor_code(record: &Ribonucleotide) -> Result<String> {
    // Source stores the first raw-code byte, not its origin or display text.
    // A non-ASCII leading byte is not a valid standalone Rust UTF-8 string.
    match record.code().as_bytes().first() {
        Some(byte) if byte.is_ascii() => Ok(char::from(*byte).to_string()),
        _ => Err(Error::Unsupported(
            "RNA graph neighbor requires an ASCII leading code byte".into(),
        )),
    }
}

type Configuration = (
    Option<Arc<Ribonucleotide>>,
    Option<Arc<Ribonucleotide>>,
    Vec<Pattern>,
    Vec<Pattern>,
);
fn configure(enzyme: &DigestionEnzymeRNA, registry: &RibonucleotideDB) -> Result<Configuration> {
    // The inherited main regex is compiled by C++ even though RNA digest never
    // uses it. Support the same native subset here rather than silently accepting it.
    Pattern::compile(enzyme.regex())?;
    let five = match enzyme.five_prime_gain() {
        "" => None,
        "p" => Some(registry.get("5'-p")?),
        code => Some(registry.get(code)?),
    };
    let three = match enzyme.three_prime_gain() {
        "" => None,
        "p" => Some(registry.get("3'-p")?),
        "c" => Some(registry.get("3'-c")?),
        code => Some(registry.get(&format!("[{code}]"))?),
    };
    Ok((
        five,
        three,
        compile_patterns(enzyme.cuts_after())?,
        compile_patterns(enzyme.cuts_before())?,
    ))
}
struct Budget(usize);
impl Budget {
    fn tree_key(&mut self, key: &str, entries: usize) -> Result<()> {
        // Binary depth bounds B-tree height; 24 comparisons/moves per level
        // conservatively covers the current stdlib's 11-key nodes and rebalancing.
        let depth = entries.saturating_add(1).ilog2() as usize + 1;
        let work = key
            .len()
            .checked_add(1)
            .and_then(|n| n.checked_mul(depth))
            .and_then(|n| n.checked_mul(24))
            .ok_or_else(|| invalid("enzyme index work overflow"))?;
        self.charge(work)
    }
    fn charge(&mut self, amount: usize) -> Result<()> {
        self.0 = self
            .0
            .checked_sub(amount)
            .ok_or_else(|| invalid("RNA digestion work limit exceeded"))?;
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(|| invalid("RNA size overflow"))
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
