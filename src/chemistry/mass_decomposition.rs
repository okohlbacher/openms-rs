// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Amino-acid count values from the source MassDecomposition class.
//!
//! This is a composition container, not the separate mass-decomposition solver.
//! The source's cached historical maximum and map-only comparison are retained;
//! see `docs/MASS_DECOMPOSITION_SUPPORT.md` for parsing and ordering details.

use crate::data_structures::list::ListParse;
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::str::FromStr;

pub const MAX_MASS_DECOMPOSITION_INPUT_BYTES: usize = 1024 * 1024;
pub const MAX_MASS_DECOMPOSITION_OUTPUT_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_MASS_DECOMPOSITION_WORK: usize = 50_000_000;

/// Source count map plus its historical maximum. Equal maps can have different
/// cached maxima after duplicate input symbols, so [`Self::source_cmp`] is
/// deliberately separate from full-value Rust equality. No `Ord` is implemented.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MassDecomposition {
    counts: BTreeMap<u8, usize>,
    maximum: usize,
}

fn invalid(message: &'static str) -> Error {
    Error::InvalidValue(message.into())
}
fn charge(work: &mut usize, amount: usize) -> Result<()> {
    *work = work
        .checked_sub(amount)
        .ok_or_else(|| invalid("mass-decomposition work limit exceeded"))?;
    Ok(())
}
fn input(text: &str, work: &mut usize) -> Result<()> {
    if text.len() > MAX_MASS_DECOMPOSITION_INPUT_BYTES {
        return Err(invalid("mass-decomposition input exceeds one MiB"));
    }
    charge(work, text.len())
}
fn trim(text: &str) -> &str {
    text.trim_matches([' ', '\t', '\r', '\n'])
}

impl MassDecomposition {
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse literal-space-separated byte-symbol/i32-count tokens. A suffix
    /// beginning at the first '(' is ignored, and only that case trims input.
    /// Duplicate symbols keep the last count but the greatest observed maximum.
    /// Negative counts are checked errors instead of unsigned Size conversion.
    pub fn parse(text: &str) -> Result<Self> {
        let mut work = MAX_MASS_DECOMPOSITION_WORK;
        Self::parse_with_work(text, &mut work)
    }
    fn parse_with_work(text: &str, work: &mut usize) -> Result<Self> {
        input(text, work)?;
        let text = match text.find('(') {
            Some(index) => trim(&text[..index]),
            None => text,
        };
        charge(work, text.len())?;
        let mut result = Self::new();
        if text.is_empty() {
            return Ok(result);
        }
        for token in text.split(' ') {
            // At most128 byte keys exist: this conservatively charges every
            // possible key comparison before each BTreeMap insertion.
            charge(work, result.counts.len() + token.len() + 1)?;
            let symbol = *token
                .as_bytes()
                .first()
                .ok_or_else(|| invalid("empty mass-decomposition token"))?;
            if !symbol.is_ascii() {
                return Err(invalid(
                    "mass-decomposition symbols must be single ASCII bytes",
                ));
            }
            let count = i32::from_list_item(&token[1..])?;
            let count =
                usize::try_from(count).map_err(|_| invalid("negative mass-decomposition count"))?;
            result.maximum = result.maximum.max(count);
            result.counts.insert(symbol, count);
        }
        Ok(result)
    }

    /// Compact sorted map, with explicit zero counts and source outer trimming.
    /// This text may not reconstruct a historical maximum from duplicate input.
    pub fn to_text(&self) -> Result<String> {
        let mut length = 0usize;
        for count in self.counts.values() {
            let digits = if *count == 0 {
                1
            } else {
                count.ilog10() as usize + 1
            };
            length = length
                .checked_add(digits + 2)
                .ok_or_else(|| invalid("decomposition text size overflow"))?;
        }
        if length > MAX_MASS_DECOMPOSITION_OUTPUT_BYTES {
            return Err(invalid("decomposition output byte limit exceeded"));
        }
        let mut text = String::with_capacity(length);
        use std::fmt::Write;
        for (&symbol, count) in &self.counts {
            write!(text, "{}{} ", char::from(symbol), count)
                .expect("writing to String cannot fail");
        }
        let begin = text.len() - text.trim_start_matches([' ', '\t', '\r', '\n']).len();
        let end = text.trim_end_matches([' ', '\t', '\r', '\n']).len();
        text.truncate(end);
        text.drain(..begin);
        Ok(text)
    }

    /// Repeat symbols in ascending byte order, checking total bytes before any
    /// count-sized allocation or loop. Zero-count keys produce no characters.
    pub fn to_expanded_string(&self) -> Result<String> {
        let mut length = 0usize;
        for &count in self.counts.values() {
            length = length
                .checked_add(count)
                .ok_or_else(|| invalid("expanded decomposition size overflow"))?;
        }
        if length > MAX_MASS_DECOMPOSITION_OUTPUT_BYTES {
            return Err(invalid("expanded decomposition exceeds16MiB"));
        }
        let mut work = MAX_MASS_DECOMPOSITION_WORK;
        charge(&mut work, self.counts.len() + length)?;
        let mut bytes = Vec::with_capacity(length);
        for (&symbol, &count) in &self.counts {
            bytes.resize(bytes.len() + count, symbol);
        }
        Ok(String::from_utf8(bytes).expect("decomposition keys are ASCII"))
    }

    /// Source operator+: new keys compare their counts with the original left
    /// maximum, so a later new key can lower the result's cached maximum.
    /// Existing keys compare their sums with the current result maximum.
    pub fn checked_add(&self, rhs: &Self) -> Result<Self> {
        let mut result = self.clone();
        for (&symbol, &count) in &rhs.counts {
            if let Some(existing) = result.counts.get_mut(&symbol) {
                *existing = existing
                    .checked_add(count)
                    .ok_or_else(|| invalid("mass-decomposition count overflow"))?;
                result.maximum = result.maximum.max(*existing);
            } else {
                result.counts.insert(symbol, count);
                if count > self.maximum {
                    result.maximum = count;
                }
            }
        }
        Ok(result)
    }

    /// Add counts without changing self if any unsigned count would overflow.
    pub fn checked_add_assign(&mut self, rhs: &Self) -> Result<()> {
        // Counts are checked before committing any partial addition.
        let mut maximum = self.maximum;
        for (&symbol, &count) in &rhs.counts {
            let sum = self
                .counts
                .get(&symbol)
                .copied()
                .unwrap_or(0)
                .checked_add(count)
                .ok_or_else(|| invalid("mass-decomposition count overflow"))?;
            maximum = maximum.max(sum);
        }
        for (&symbol, &count) in &rhs.counts {
            *self.counts.entry(symbol).or_default() += count;
        }
        self.maximum = maximum;
        Ok(())
    }

    /// Cached maximum, including earlier counts overwritten by duplicate tokens.
    pub fn number_of_max_aa(&self) -> usize {
        self.maximum
    }

    /// Source map ordering: lexicographic symbol/count pairs, ignoring maximum.
    pub fn source_cmp(&self, rhs: &Self) -> Ordering {
        self.counts.cmp(&rhs.counts)
    }
    pub fn source_less(&self, rhs: &Self) -> bool {
        self.source_cmp(rhs).is_lt()
    }

    /// Source equality against parsed text compares both map and cached maximum.
    pub fn equals_text(&self, text: &str) -> Result<bool> {
        Ok(self == &Self::parse(text)?)
    }

    /// Multiset containment: tag order is ignored; every byte needs enough count.
    /// An absent symbol returns false, including non-ASCII UTF-8 tag bytes.
    pub fn contains_tag(&self, tag: &str) -> Result<bool> {
        let mut work = MAX_MASS_DECOMPOSITION_WORK;
        input(tag, &mut work)?;
        let mut used = BTreeMap::<u8, usize>::new();
        for symbol in tag.bytes() {
            charge(&mut work, self.counts.len() + used.len() + 1)?;
            let Some(&available) = self.counts.get(&symbol) else {
                return Ok(false);
            };
            let count = used.entry(symbol).or_default();
            *count += 1; // Tag length is bounded well below usize::MAX.
            if *count > available {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// True when every key/count in rhs exists here with at least that count.
    /// A missing key with count zero still fails, unlike an empty tag.
    /// No source stderr diagnostic is emitted on an incompatible composition.
    pub fn compatible(&self, rhs: &Self) -> bool {
        rhs.counts.iter().all(|(symbol, count)| {
            self.counts
                .get(symbol)
                .is_some_and(|available| available >= count)
        })
    }
}
impl FromStr for MassDecomposition {
    type Err = Error;
    fn from_str(text: &str) -> Result<Self> {
        Self::parse(text)
    }
}
