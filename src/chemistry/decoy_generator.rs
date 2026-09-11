// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Seeded protein and peptide decoys from pinned OpenMS `7c029e8`.
//!
//! Modifications are rejected. Peptide anchors, repeated in-place shuffles and
//! the text-only cache follow the source, including its positional exceptions.
//! Unspecific cleavage concatenates overlapping products, so those outputs are
//! not isobaric with the input. Random operations have checked cumulative limits.

use super::{AASequence, Protease, SequenceModification};
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::mem::size_of;
use std::ops::Range;
use std::time::{SystemTime, UNIX_EPOCH};

#[path = "decoy_random.rs"]
mod decoy_random;
pub(crate) use decoy_random::DecoyRandom;

pub const MAX_DECOY_INPUT_RESIDUES: usize = 1_000_000;
pub const MAX_DECOY_OUTPUT_RESIDUES: usize = 1_000_000;
pub const MAX_DECOY_PEPTIDES: usize = 100_000;
pub const MAX_DECOY_VARIANTS: usize = 1_000;
pub const MAX_DECOY_WORK: usize = 50_000_000;
pub const MAX_DECOY_BYTES: usize = 256 * 1024 * 1024;
pub const MAX_DECOY_CACHE_ENTRIES: usize = 100_000;
pub const MAX_DECOY_CACHE_BYTES: usize = 16 * 1024 * 1024;
// Includes a conservatively sparse B-tree node, not just the two String headers.
const CACHE_ENTRY_BYTES: usize = 1024;

/// Stateful peptide shuffler. Cache keys contain only the unmodified peptide
/// text: enzyme, positional anchor, attempt count and seed are deliberately not
/// part of the key. Reseeding therefore leaves existing cached choices intact.
///
/// `shuffle_peptides` stages its RNG and new cache entries, committing only on
/// success. `shuffle` instead creates the source's independent per-product RNGs
/// and never changes this receiver. Concurrent mutation requires caller locking.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecoyGenerator {
    random: DecoyRandom,
    cache: BTreeMap<String, String>,
    cache_bytes: usize,
}

impl Default for DecoyGenerator {
    fn default() -> Self {
        let ticks = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_else(|error| error.duration())
            .as_nanos();
        Self::with_seed(ticks as u64)
    }
}

impl DecoyGenerator {
    /// Reproducible constructor; unlike `Default`, does not read the clock.
    pub fn with_seed(seed: u64) -> Self {
        Self {
            random: DecoyRandom::seeded(seed),
            cache: BTreeMap::new(),
            cache_bytes: 0,
        }
    }

    /// Reset the random stream without clearing the source peptide cache.
    pub fn set_seed(&mut self, seed: u64) {
        self.random.reseed(seed);
    }

    /// Reverse every parent residue. Empty input stays empty; no state changes.
    pub fn reverse_protein(&self, protein: &AASequence) -> Result<AASequence> {
        let mut work = DecoyWork::default();
        validate_protein(protein, &mut work)?;
        work.consume(protein.len())?;
        work.allocate(protein.len())?;
        let mut bytes = protein.as_str().as_bytes().to_vec();
        bytes.reverse();
        finish_sequence(&bytes, &mut work)
    }

    /// Reverse each digested product except its final residue, then fully reverse
    /// the final product. This source positional rule also applies to N-terminal
    /// enzymes and proteins ending at a cleavage site. Empty input is an error.
    pub fn reverse_peptides(&self, protein: &AASequence, enzyme: Protease) -> Result<AASequence> {
        let mut work = DecoyWork::default();
        validate_protein(protein, &mut work)?;
        let ranges = digest_ranges(protein.as_str(), enzyme, &mut work)?;
        require_products(&ranges)?;
        let length = products_length(&ranges, &mut work)?;
        work.check_output(length)?;
        work.allocate(length)?;
        work.consume(length)?;
        let mut output = Vec::with_capacity(length);
        for (index, range) in ranges.iter().enumerate() {
            let peptide = &protein.as_str().as_bytes()[range.clone()];
            if index + 1 == ranges.len() {
                output.extend(peptide.iter().rev());
            } else {
                output.extend(peptide[..peptide.len() - 1].iter().rev());
                output.push(peptide[peptide.len() - 1]);
            }
        }
        finish_sequence(&output, &mut work)
    }

    /// Minimize the maximum forward/reverse positional identity. Each attempt
    /// permutes the previous candidate in place; exact ties keep the earlier
    /// best. Zero attempts cache the original products. Empty input is an error.
    /// Source negative attempt counts have the same behavior as native zero.
    pub fn shuffle_peptides(
        &mut self,
        protein: &AASequence,
        enzyme: Protease,
        max_attempts: usize,
    ) -> Result<AASequence> {
        self.shuffle_peptides_with_work(protein, enzyme, max_attempts, &mut DecoyWork::default())
    }

    fn shuffle_peptides_with_work(
        &mut self,
        protein: &AASequence,
        enzyme: Protease,
        max_attempts: usize,
        work: &mut DecoyWork,
    ) -> Result<AASequence> {
        validate_protein(protein, work)?;
        work.consume(size_of::<DecoyRandom>() / size_of::<u64>())?;
        work.allocate(size_of::<DecoyRandom>())?;
        let mut random = self.random.clone();
        let mut delta = CacheDelta::default();
        let output = shuffle_raw(
            protein.as_str(),
            enzyme,
            max_attempts,
            &self.cache,
            self.cache_bytes,
            &mut delta,
            &mut random,
            work,
        )?;
        let result = finish_sequence(&output, work)?;
        // Cache commit comparisons were precharged when each delta entry was
        // created. Move only new entries, never clone or traverse the old cache.
        self.cache.extend(delta.entries);
        self.cache_bytes += delta.bytes;
        self.random = random;
        Ok(result)
    }

    /// Build independent complete variants. Every outer product longer than two
    /// residues gets a fresh local shuffler seeded with `4711 + variant`, which
    /// redigests that product. Receiver RNG/cache are unused. Source default
    /// factor is one; native zero represents source nonpositive factors.
    /// Positive factors with empty input produce that many empty sequences.
    pub fn shuffle(
        &self,
        protein: &AASequence,
        enzyme: Protease,
        decoy_factor: usize,
    ) -> Result<Vec<AASequence>> {
        let mut work = DecoyWork::default();
        validate_protein(protein, &mut work)?;
        let ranges = digest_ranges(protein.as_str(), enzyme, &mut work)?;
        if decoy_factor > MAX_DECOY_VARIANTS {
            return Err(invalid("decoy variant limit exceeded"));
        }
        // Unspecific inner digestion expands each longer outer product again.
        // Compute that complete output size before any candidate/RNG work.
        let mut variant_length = 0;
        work.consume(ranges.len())?;
        for range in &ranges {
            let length = if enzyme == Protease::UnspecificCleavage && range.len() > 2 {
                unspecific_residues(range.len())?
            } else {
                range.len()
            };
            variant_length = add(variant_length, length)?;
        }
        work.check_output(multiply(variant_length, decoy_factor)?)?;
        work.allocate(multiply(decoy_factor, size_of::<AASequence>())?)?;
        let mut results = Vec::with_capacity(decoy_factor);
        for variant in 0..decoy_factor {
            work.consume(1)?;
            work.allocate(variant_length)?;
            let mut output = Vec::with_capacity(variant_length);
            for range in &ranges {
                let peptide = &protein.as_str()[range.clone()];
                if peptide.len() <= 2 {
                    work.consume(peptide.len())?;
                    output.extend_from_slice(peptide.as_bytes());
                } else {
                    // Source intentionally resets seed/cache for every product.
                    work.consume(size_of::<DecoyRandom>() / size_of::<u64>())?;
                    work.allocate(size_of::<DecoyRandom>())?;
                    let mut random = DecoyRandom::seeded(4711 + variant as u64);
                    let mut delta = CacheDelta::default();
                    let shuffled = shuffle_raw(
                        peptide,
                        enzyme,
                        100,
                        &BTreeMap::new(),
                        0,
                        &mut delta,
                        &mut random,
                        &mut work,
                    )?;
                    work.consume(shuffled.len())?;
                    output.extend_from_slice(&shuffled);
                }
            }
            results.push(finish_sequence(&output, &mut work)?);
        }
        Ok(results)
    }
}

#[derive(Default)]
struct CacheDelta {
    entries: BTreeMap<String, String>,
    bytes: usize,
}

// A single budget follows all inner digestions, attempts, rejected random draws,
// cache comparisons/copies and output construction in one public operation.
struct DecoyWork {
    remaining: usize,
    bytes: usize,
    peptides: usize,
    output: usize,
}
impl Default for DecoyWork {
    fn default() -> Self {
        Self {
            remaining: MAX_DECOY_WORK,
            bytes: MAX_DECOY_BYTES,
            peptides: MAX_DECOY_PEPTIDES,
            output: MAX_DECOY_OUTPUT_RESIDUES,
        }
    }
}
impl DecoyWork {
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("decoy cumulative work limit exceeded"))?;
        Ok(())
    }
    fn allocate(&mut self, bytes: usize) -> Result<()> {
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("decoy cumulative allocation limit exceeded"))?;
        Ok(())
    }
    fn check_output(&self, count: usize) -> Result<()> {
        if count > self.output {
            return Err(invalid("decoy cumulative output residue limit exceeded"));
        }
        Ok(())
    }
    fn products(&mut self, count: usize) -> Result<()> {
        self.peptides = self
            .peptides
            .checked_sub(count)
            .ok_or_else(|| invalid("decoy cumulative peptide product limit exceeded"))?;
        self.consume(count)
    }
    fn cache_search(&mut self, entries: usize, length: usize) -> Result<()> {
        // B-tree nodes perform several lexicographic comparisons per level.
        // Charge a conservative 32 comparisons per binary-tree level, covering
        // the delta/base lookup and a later insertion without old-cache scans.
        let levels = usize::BITS as usize - entries.leading_zeros() as usize + 1;
        self.consume(multiply(multiply(32, levels)?, add(length, 1)?)?)
    }
}

fn validate_protein(protein: &AASequence, work: &mut DecoyWork) -> Result<()> {
    if protein.len() > MAX_DECOY_INPUT_RESIDUES {
        return Err(invalid("decoy input residue limit exceeded"));
    }
    work.consume(add(protein.len(), 1)?)?;
    if protein.is_modified() {
        return Err(Error::Unsupported(
            "decoy generation requires an unmodified protein".into(),
        ));
    }
    Ok(())
}

fn digest_ranges(
    sequence: &str,
    enzyme: Protease,
    work: &mut DecoyWork,
) -> Result<Vec<Range<usize>>> {
    let n = sequence.len();
    if enzyme == Protease::UnspecificCleavage {
        let count = usize::try_from((n as u128) * (n as u128 + 1) / 2)
            .map_err(|_| invalid("decoy digestion count overflow"))?;
        work.products(count)?;
        work.allocate(multiply(count, size_of::<Range<usize>>())?)?;
        let mut ranges = Vec::with_capacity(count);
        // AASequence digestion overrides missed cleavages for this enzyme.
        // Its source order is length, then start, unlike the native generic API.
        for length in 1..=n {
            for start in 0..=n - length {
                ranges.push(start..start + length);
            }
        }
        Ok(ranges)
    } else {
        work.consume(add(multiply(n, 3)?, 1)?)?;
        // cleavage_sites validates/scans once and grows a vector; four times its
        // maximum payload bounds geometric growth plus the live final buffer.
        work.allocate(multiply(add(n, 1)?, 4 * size_of::<usize>())?)?;
        let cuts = enzyme.cleavage_sites(sequence)?;
        let count = cuts.len().saturating_sub(1);
        work.products(count)?;
        work.allocate(multiply(count, size_of::<Range<usize>>())?)?;
        Ok(cuts.windows(2).map(|cut| cut[0]..cut[1]).collect())
    }
}

fn require_products(ranges: &[Range<usize>]) -> Result<()> {
    if ranges.is_empty() {
        return Err(invalid(
            "peptide decoy generation requires a nonempty protein",
        ));
    }
    Ok(())
}
fn products_length(ranges: &[Range<usize>], work: &mut DecoyWork) -> Result<usize> {
    work.consume(ranges.len())?;
    ranges
        .iter()
        .try_fold(0, |sum, range| add(sum, range.len()))
}
fn unspecific_residues(n: usize) -> Result<usize> {
    usize::try_from((n as u128) * (n as u128 + 1) * (n as u128 + 2) / 6)
        .map_err(|_| invalid("unspecific decoy output length overflow"))
}

#[allow(clippy::too_many_arguments)]
fn shuffle_raw(
    sequence: &str,
    enzyme: Protease,
    max_attempts: usize,
    cache: &BTreeMap<String, String>,
    cache_bytes: usize,
    delta: &mut CacheDelta,
    random: &mut DecoyRandom,
    work: &mut DecoyWork,
) -> Result<Vec<u8>> {
    let ranges = digest_ranges(sequence, enzyme, work)?;
    require_products(&ranges)?;
    let length = products_length(&ranges, work)?;
    work.check_output(length)?;
    work.allocate(length)?;
    let mut output = Vec::with_capacity(length);
    for (index, range) in ranges.iter().enumerate() {
        let target = &sequence[range.clone()];
        let entries = add(cache.len(), delta.entries.len())?;
        work.cache_search(entries, target.len())?;
        if let Some(cached) = delta.entries.get(target).or_else(|| cache.get(target)) {
            work.consume(cached.len())?;
            output.extend_from_slice(cached.as_bytes());
            continue;
        }
        let entry_bytes = add(multiply(target.len(), 2)?, CACHE_ENTRY_BYTES)?;
        if entries >= MAX_DECOY_CACHE_ENTRIES
            || add(add(cache_bytes, delta.bytes)?, entry_bytes)? > MAX_DECOY_CACHE_BYTES
        {
            return Err(invalid("decoy persistent peptide cache limit exceeded"));
        }
        // Candidate + best buffers, owned key/value, and conservatively sparse
        // nodes in both the delta and the eventual destination cache.
        work.allocate(add(
            add(multiply(target.len(), 2)?, entry_bytes)?,
            CACHE_ENTRY_BYTES,
        )?)?;
        work.consume(multiply(target.len(), 2)?)?;
        let mut candidate = target.as_bytes().to_vec();
        let mut best = candidate.clone();
        let mut lowest = 1.0;
        let final_product = index + 1 == ranges.len();
        let shuffle_end = target.len() - usize::from(!final_product);
        for _ in 0..max_attempts {
            work.consume(add(multiply(target.len(), 2)?, 1)?)?;
            random.shuffle(&mut candidate[..shuffle_end], &mut work.remaining)?;
            let identity = sequence_identity(&candidate, target.as_bytes());
            if identity < lowest {
                lowest = identity;
                work.consume(target.len())?;
                best.copy_from_slice(&candidate);
                if (final_product && identity == 0.0)
                    || (!final_product && identity <= 1.0 / target.len() as f64 + 1e-6)
                {
                    break;
                }
            }
        }
        work.consume(multiply(target.len(), 3)?)?;
        work.cache_search(add(entries, 1)?, target.len())?;
        output.extend_from_slice(&best);
        // All source letters are ASCII and permutation cannot introduce syntax.
        let value = String::from_utf8(best)
            .map_err(|_| invalid("internal decoy sequence encoding error"))?;
        delta.entries.insert(target.to_owned(), value);
        delta.bytes += entry_bytes;
    }
    Ok(output)
}

fn sequence_identity(decoy: &[u8], target: &[u8]) -> f64 {
    debug_assert_eq!(decoy.len(), target.len());
    debug_assert!(!target.is_empty());
    let forward = target.iter().zip(decoy).filter(|(a, b)| a == b).count();
    let reverse = target
        .iter()
        .zip(decoy.iter().rev())
        .filter(|(a, b)| a == b)
        .count();
    (forward as f64 / target.len() as f64).max(reverse as f64 / target.len() as f64)
}

fn finish_sequence(bytes: &[u8], work: &mut DecoyWork) -> Result<AASequence> {
    work.check_output(bytes.len())?;
    // The unmodified parser fills a parent string and optional annotation slots,
    // with a bounded small-element formula cache. Charge geometric buffer growth
    // and the parser/chemistry residue scans before constructing the result.
    work.consume(add(multiply(bytes.len(), 16)?, 1)?)?;
    work.allocate(add(
        multiply(
            bytes.len(),
            4 * (size_of::<Option<SequenceModification>>() + 1),
        )?,
        4096,
    )?)?;
    let text = std::str::from_utf8(bytes)
        .map_err(|_| invalid("internal decoy sequence encoding error"))?;
    let sequence = AASequence::parse(text)?;
    work.output -= bytes.len();
    Ok(sequence)
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}
fn add(left: usize, right: usize) -> Result<usize> {
    left.checked_add(right)
        .ok_or_else(|| invalid("decoy resource accounting overflow"))
}
fn multiply(left: usize, right: usize) -> Result<usize> {
    left.checked_mul(right)
        .ok_or_else(|| invalid("decoy resource accounting overflow"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_counts_forward_and_reverse_and_keeps_the_larger_fraction() {
        assert_eq!(sequence_identity(b"ABCD", b"ABCD"), 1.0);
        assert_eq!(sequence_identity(b"DCBA", b"ABCD"), 1.0);
        assert_eq!(sequence_identity(b"BADC", b"ABCD"), 0.0);
        assert_eq!(sequence_identity(b"ACBD", b"ABCD"), 0.5);
        assert_eq!(sequence_identity(b"AAAA", b"AAAA"), 1.0);
    }

    #[test]
    fn late_work_exhaustion_after_a_cache_miss_does_not_commit_rng_or_delta() {
        let mut generator = DecoyGenerator::with_seed(4711);
        generator
            .shuffle_peptides(
                &AASequence::parse("PEPTIDE").unwrap(),
                Protease::NoCleavage,
                3,
            )
            .unwrap();
        let original = generator.clone();
        let protein = AASequence::parse(&format!("ACDR{}", "A".repeat(50))).unwrap();
        let mut work = DecoyWork {
            remaining: 5000,
            ..Default::default()
        };
        let error = generator
            .shuffle_peptides_with_work(&protein, Protease::TrypsinP, 100, &mut work)
            .unwrap_err();
        assert!(error.to_string().contains("work limit"));
        assert_eq!(generator, original);
    }

    #[test]
    fn cache_limit_and_native_output_allocation_failure_are_transactional() {
        let mut generator = DecoyGenerator::with_seed(8);
        generator.cache_bytes = MAX_DECOY_CACHE_BYTES;
        let original = generator.clone();
        let protein = AASequence::parse("ACDE").unwrap();
        assert!(
            generator
                .shuffle_peptides(&protein, Protease::NoCleavage, 1)
                .is_err()
        );
        assert_eq!(generator, original);

        let mut generator = DecoyGenerator::with_seed(8);
        let original = generator.clone();
        let mut work = DecoyWork {
            bytes: 6500,
            ..Default::default()
        };
        assert!(
            generator
                .shuffle_peptides_with_work(&protein, Protease::NoCleavage, 0, &mut work)
                .unwrap_err()
                .to_string()
                .contains("allocation limit")
        );
        assert_eq!(generator, original);
    }
}
