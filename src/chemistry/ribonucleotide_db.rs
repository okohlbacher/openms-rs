// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Ordered RNA providers and immutable registries from OpenMS4-core 7c029e8.
//! Iteration retains duplicates; exact lookup and ambiguity targets use the last
//! record for a code. Unresolved later ambiguities retain earlier valid mappings.

use super::ribonucleotide::{
    MAX_RIBONUCLEOTIDE_CODE_BYTES, Ribonucleotide, RibonucleotideRecord,
    RibonucleotideTermSpecificity,
};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::{Arc, OnceLock};

#[path = "ribonucleotide_providers.rs"]
mod providers;
#[cfg(feature = "rna-json")]
pub use providers::read_modomics_json;
pub use providers::read_tsv;

pub(super) const MAX_RECORDS: usize = 100_000;
pub(super) const MAX_BYTES: usize = 128 * 1024 * 1024;
pub(super) const MAX_WORK: usize = 50_000_000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RibonucleotideEntry {
    pub ribonucleotide: Arc<Ribonucleotide>,
    pub alternatives: Option<[String; 2]>,
}
impl RibonucleotideEntry {
    /// Source entry predicate depends on the first alternative, not code suffix.
    pub fn is_ambiguous(&self) -> bool {
        self.alternatives.as_ref().is_some_and(|a| !a[0].is_empty())
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RibonucleotideDiagnostic {
    /// One-based TSV line, JSON iteration index, or registry input-entry index.
    pub index: usize,
    pub message: String,
}

#[derive(Clone, Debug, Default)]
pub struct RibonucleotideLoadReport {
    pub entries: Vec<RibonucleotideEntry>,
    pub skipped: usize,
    pub diagnostics: Vec<RibonucleotideDiagnostic>,
}

/// Owned immutable registry. Construction allows 100,000 records, 128 MiB of
/// conservative logical storage, and 50M charged key/work units. Handles survive
/// registry destruction; cloning the registry shares its immutable records.
#[derive(Clone, Debug, Default)]
pub struct RibonucleotideDB {
    entries: Vec<Arc<Ribonucleotide>>,
    by_code: BTreeMap<String, usize>,
    code_lengths: BTreeSet<usize>,
    ambiguity: BTreeMap<String, [usize; 2]>,
    diagnostics: Vec<RibonucleotideDiagnostic>,
}

impl RibonucleotideDB {
    pub fn from_records(records: Vec<Ribonucleotide>) -> Result<Self> {
        if records.len() > MAX_RECORDS {
            return Err(invalid("record limit exceeded"));
        }
        let mut bytes = 0;
        for record in &records {
            add(
                &mut bytes,
                record.payload_bytes()?,
                MAX_BYTES,
                "registry bytes",
            )?;
        }
        Self::from_entries(
            records
                .into_iter()
                .map(|ribonucleotide| RibonucleotideEntry {
                    ribonucleotide: Arc::new(ribonucleotide),
                    alternatives: None,
                })
                .collect(),
        )
    }

    pub fn from_entries(entries: Vec<RibonucleotideEntry>) -> Result<Self> {
        if entries.len() > MAX_RECORDS {
            return Err(invalid("record limit exceeded"));
        }
        let mut bytes = entries
            .len()
            .saturating_mul(std::mem::size_of::<RibonucleotideEntry>());
        let mut work = 0;
        for entry in &entries {
            add(
                &mut bytes,
                entry.ribonucleotide.payload_bytes()?,
                MAX_BYTES,
                "registry bytes",
            )?;
            add(
                &mut bytes,
                entry.ribonucleotide.code().len() + 128,
                MAX_BYTES,
                "registry bytes",
            )?;
            key_work(&mut work, entry.ribonucleotide.code())?;
            if let Some(alternatives) = &entry.alternatives {
                for alternative in alternatives {
                    if alternative.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES {
                        return Err(invalid("alternative code byte limit exceeded"));
                    }
                    add(
                        &mut bytes,
                        alternative.len() + 64,
                        MAX_BYTES,
                        "registry bytes",
                    )?;
                    key_work(&mut work, alternative)?;
                }
            }
        }
        let mut result = Self::default();
        result
            .entries
            .try_reserve_exact(entries.len())
            .map_err(|_| invalid("registry allocation failed"))?;
        for entry in &entries {
            result
                .by_code
                .insert(entry.ribonucleotide.code().into(), result.entries.len());
            result
                .code_lengths
                .insert(entry.ribonucleotide.code().len());
            result.entries.push(Arc::clone(&entry.ribonucleotide));
        }
        // Resolve after every provider has contributed to the final winner map.
        for (index, entry) in entries.iter().enumerate() {
            if let Some(alternatives) = entry.alternatives.as_ref().filter(|a| !a[0].is_empty()) {
                let first = result.by_code.get(&alternatives[0]);
                let second = result.by_code.get(&alternatives[1]);
                let code = entry.ribonucleotide.code();
                key_work(&mut work, code)?;
                if let (Some(&first), Some(&second)) = (first, second) {
                    add(&mut bytes, code.len() + 96, MAX_BYTES, "registry bytes")?;
                    result.ambiguity.insert(code.into(), [first, second]);
                } else {
                    let message = format!("unresolved alternatives for ribonucleotide {code:?}");
                    add(
                        &mut bytes,
                        message.len() + 64,
                        MAX_BYTES,
                        "registry diagnostics",
                    )?;
                    result.diagnostics.push(RibonucleotideDiagnostic {
                        index: index + 1,
                        message,
                    });
                }
            }
        }
        Ok(result)
    }

    pub fn entries(&self) -> &[Arc<Ribonucleotide>] {
        &self.entries
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn diagnostics(&self) -> &[RibonucleotideDiagnostic] {
        &self.diagnostics
    }

    pub fn get(&self, code: &str) -> Result<Arc<Ribonucleotide>> {
        if code.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES {
            return Err(invalid("lookup code byte limit exceeded"));
        }
        self.by_code
            .get(code)
            .map(|&index| Arc::clone(&self.entries[index]))
            .ok_or_else(|| invalid("ribonucleotide code was not found"))
    }

    /// Longest complete UTF-8 code prefix, independent of sequence tokenization.
    pub fn get_prefix(&self, sequence: &str) -> Result<Arc<Ribonucleotide>> {
        let mut work = 0;
        for &length in self.code_lengths.range(..=sequence.len()).rev() {
            if let Some(prefix) = sequence.get(..length) {
                key_work(&mut work, prefix)?;
                if let Some(&index) = self.by_code.get(prefix) {
                    return Ok(Arc::clone(&self.entries[index]));
                }
            }
        }
        Err(invalid("ribonucleotide prefix was not found"))
    }

    pub fn alternatives(&self, code: &str) -> Result<[Arc<Ribonucleotide>; 2]> {
        if code.len() > MAX_RIBONUCLEOTIDE_CODE_BYTES {
            return Err(invalid("lookup code byte limit exceeded"));
        }
        self.ambiguity
            .get(code)
            .map(|&[a, b]| [Arc::clone(&self.entries[a]), Arc::clone(&self.entries[b])])
            .ok_or_else(|| invalid("ribonucleotide alternatives were not found"))
    }

    pub fn global() -> &'static Self {
        static DB: OnceLock<RibonucleotideDB> = OnceLock::new();
        DB.get_or_init(|| {
            let records: &[EmbeddedRibonucleotide] =
                include!("../../resources/rna/ribonucleotides.rs");
            let entries = records
                .iter()
                .map(|record| {
                    let ribonucleotide = Ribonucleotide::from_record(RibonucleotideRecord {
                        name: record.name.into(),
                        code: record.code.into(),
                        new_code: record.new_code.into(),
                        html_code: record.html_code.into(),
                        formula: record.formula.parse()?,
                        origin: record.origin,
                        mono_mass: f64::from_bits(record.mono_mass_bits),
                        average_mass: f64::from_bits(record.average_mass_bits),
                        term_specificity: match record.term_specificity {
                            0 => RibonucleotideTermSpecificity::Anywhere,
                            1 => RibonucleotideTermSpecificity::FivePrime,
                            2 => RibonucleotideTermSpecificity::ThreePrime,
                            _ => return Err(invalid("invalid embedded terminal specificity")),
                        },
                        baseloss_formula: record.baseloss_formula.parse()?,
                    })?;
                    Ok(RibonucleotideEntry {
                        ribonucleotide: Arc::new(ribonucleotide),
                        alternatives: record.alternatives.map(|[a, b]| [a.into(), b.into()]),
                    })
                })
                .collect::<Result<Vec<_>>>()
                .expect("validated embedded RNA records");
            Self::from_entries(entries).expect("validated embedded RNA registry")
        })
    }
}

struct EmbeddedRibonucleotide {
    name: &'static str,
    code: &'static str,
    new_code: &'static str,
    html_code: &'static str,
    formula: &'static str,
    origin: char,
    mono_mass_bits: u64,
    average_mass_bits: u64,
    term_specificity: u8,
    baseloss_formula: &'static str,
    alternatives: Option<[&'static str; 2]>,
}

pub(super) fn add(total: &mut usize, amount: usize, limit: usize, label: &str) -> Result<()> {
    *total = total
        .checked_add(amount)
        .filter(|&n| n <= limit)
        .ok_or_else(|| invalid(&format!("{label} limit exceeded")))?;
    Ok(())
}
fn key_work(work: &mut usize, key: &str) -> Result<()> {
    // Fewer than 1024 key comparisons for a usize-sized BTree. String comparison
    // never consumes more bytes than the queried key. Charge before allocation.
    add(
        work,
        key.len().saturating_add(1).saturating_mul(1024),
        MAX_WORK,
        "registry key work",
    )
}
pub(super) fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
