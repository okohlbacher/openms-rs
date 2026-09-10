// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Bounded residue-mass sequence tags from the pinned OpenMS `Tagger`.
//! Input index order and the source's early maximum-gap pruning are retained;
//! unsorted finite coordinates are accepted without sorting them implicitly.

use super::ModificationsDB;
use crate::comparison::Tolerance;
use crate::{Error, MSSpectrum, Result};
use std::cmp::Ordering;
use std::mem::size_of;

#[path = "tagger_masses.rs"]
mod tagger_masses;

pub const MAX_TAGGER_PEAKS: usize = 1_000_000;
/// Existing entries plus emitted paths, before final deduplication.
pub const MAX_TAGGER_TAGS: usize = 1_000_000;
/// Existing string payload plus all emitted tag copies.
pub const MAX_TAGGER_TEXT_BYTES: usize = 64 * 1024 * 1024;
pub const MAX_TAGGER_WORK: usize = 50_000_000;
pub const MAX_TAGGER_BYTES: usize = 256 * 1024 * 1024;

/// All source constructor options, with explicit absolute/ppm units.
#[derive(Clone, Debug, PartialEq)]
pub struct TaggerOptions {
    pub min_tag_length: usize,
    pub max_tag_length: usize,
    pub min_charge: usize,
    pub max_charge: usize,
    pub tolerance: Tolerance,
    pub fixed_mods: Vec<String>,
    pub variable_mods: Vec<String>,
}
impl TaggerOptions {
    /// Source defaults: maximum length 65535, charge one, no modifications.
    /// Tolerance magnitude is normalized by `Tagger` construction, as in C++.
    pub fn new(min_tag_length: usize, tolerance: Tolerance) -> Self {
        Self {
            min_tag_length,
            max_tag_length: 65_535,
            min_charge: 1,
            max_charge: 1,
            tolerance,
            fixed_mods: Vec::new(),
            variable_mods: Vec::new(),
        }
    }
}

/// An owned resolved mass table. Registry lifetimes and later registry changes
/// cannot change this tagger. Modified masses still emit parent residue letters.
#[derive(Clone, Debug, PartialEq)]
pub struct Tagger {
    options: TaggerOptions,
    masses: Vec<(f64, u8)>,
    min_gap: f64,
    max_gap: f64,
}
impl Tagger {
    pub fn new(options: TaggerOptions) -> Result<Self> {
        Self::with_registry(options, ModificationsDB::global())
    }

    pub fn with_registry(mut options: TaggerOptions, registry: &ModificationsDB) -> Result<Self> {
        check_max_charge(options.max_charge)?;
        options.tolerance = match options.tolerance {
            Tolerance::Absolute(value) => Tolerance::Absolute(finite(value, "tolerance")?.abs()),
            Tolerance::Ppm(value) => Tolerance::Ppm(finite(value, "tolerance")?.abs()),
        };
        let masses = tagger_masses::resolve_mass_table(
            &options.fixed_mods,
            &options.variable_mods,
            registry,
        )?;
        let first = masses
            .first()
            .ok_or_else(|| invalid("empty tagger residue mass table"))?
            .0;
        let last = masses
            .last()
            .ok_or_else(|| invalid("empty tagger residue mass table"))?
            .0;
        let min_gap = finite(
            first - delta(options.tolerance, first)?,
            "minimum residue gap",
        )?;
        let max_gap = finite(
            last + delta(options.tolerance, last)?,
            "maximum residue gap",
        )?;
        Ok(Self {
            options,
            masses,
            min_gap,
            max_gap,
        })
    }

    /// Effective options, including the source-normalized tolerance magnitude.
    pub fn options(&self) -> &TaggerOptions {
        &self.options
    }

    /// Source setter; a value below the minimum produces an empty charge range.
    /// Zero is valid. usize::MAX is rejected because the source loop overflows.
    pub fn set_max_charge(&mut self, max_charge: usize) -> Result<()> {
        check_max_charge(max_charge)?;
        self.options.max_charge = max_charge;
        Ok(())
    }

    pub fn get_tags(&self, mzs: &[f64]) -> Result<Vec<String>> {
        let mut tags = Vec::new();
        self.append_tags(mzs, &mut tags)?;
        Ok(tags)
    }

    /// Append, lexicographically sort and deduplicate both old and new tags.
    /// The sole source early return `min_tag_length > mzs.len()` leaves the old
    /// vector entirely untouched. Any checked failure is likewise atomic.
    pub fn append_tags(&self, mzs: &[f64], tags: &mut Vec<String>) -> Result<()> {
        self.append_with_work(mzs, tags, &mut TaggerWork::default())
    }

    pub fn get_spectrum_tags(&self, spectrum: &MSSpectrum) -> Result<Vec<String>> {
        let mut tags = Vec::new();
        self.append_spectrum_tags(spectrum, &mut tags)?;
        Ok(tags)
    }

    /// Consume only m/z values. Intensities, acquisition metadata and data arrays
    /// do not affect source Tagger and are not recursively validated or copied.
    pub fn append_spectrum_tags(
        &self,
        spectrum: &MSSpectrum,
        tags: &mut Vec<String>,
    ) -> Result<()> {
        if self.options.min_tag_length > spectrum.len() {
            return Ok(());
        }
        let mut work = TaggerWork::default();
        check_peaks(spectrum.len())?;
        work.consume(spectrum.len())?;
        work.allocate::<f64>(spectrum.len())?;
        let mzs: Vec<_> = spectrum.peaks.iter().map(|peak| peak.mz).collect();
        self.append_with_work(&mzs, tags, &mut work)
    }

    fn append_with_work(
        &self,
        mzs: &[f64],
        tags: &mut Vec<String>,
        work: &mut TaggerWork,
    ) -> Result<()> {
        if self.options.min_tag_length > mzs.len() {
            return Ok(());
        }
        check_peaks(mzs.len())?;
        work.consume(mzs.len())?;
        for &mz in mzs {
            finite(mz, "input m/z")?;
        }
        work.existing_tags(tags)?;
        let mut generated = Vec::new();
        self.generate(mzs, &mut generated, work)?;
        let order = sorted_unique_indices(tags, &generated, work)?;
        work.allocate::<String>(order.len())?;
        work.consume(
            tags.len()
                .checked_add(generated.len())
                .and_then(|n| n.checked_add(order.len()))
                .ok_or_else(|| invalid("tagger commit work overflow"))?,
        )?;
        let mut result = Vec::with_capacity(order.len());
        let old_length = tags.len();
        // All failure points precede this commit. Move selected owned strings
        // instead of cloning the existing output or its potentially long labels.
        for index in order {
            result.push(if index < old_length {
                std::mem::take(&mut tags[index])
            } else {
                std::mem::take(&mut generated[index - old_length])
            });
        }
        *tags = result;
        Ok(())
    }

    fn generate(&self, mzs: &[f64], tags: &mut Vec<String>, work: &mut TaggerWork) -> Result<()> {
        let starts = mzs.len() - self.options.min_tag_length;
        if starts == 0 || self.options.min_charge > self.options.max_charge {
            return Ok(());
        }
        let depth = mzs.len().min(self.options.max_tag_length.saturating_add(1));
        work.allocate::<Frame>(depth)?;
        work.allocate::<u8>(depth)?;
        let mut stack = Vec::with_capacity(depth);
        let mut tag = Vec::with_capacity(depth);
        for start in 0..starts {
            for charge in self.options.min_charge..=self.options.max_charge {
                work.consume(1)?;
                stack.push(Frame {
                    index: start,
                    next: start + 1,
                    alternate: false,
                });
                while let Some(frame) = stack.last() {
                    work.consume(1)?;
                    let index = frame.index;
                    let next = frame.next;
                    let mut finished =
                        next == mzs.len() || tag.len() == self.options.max_tag_length;
                    let mut residue = None;
                    if !finished {
                        let gap = finite(mzs[next] - mzs[index], "peak gap")?;
                        let mass = finite(gap * charge as f64, "charged peak gap")?;
                        finished = mass > self.max_gap;
                        if !finished {
                            residue = self.amino_acid(mass, work)?;
                        }
                    }
                    if finished {
                        let mut frame = stack.pop().expect("active DFS frame");
                        if frame.alternate {
                            // Complete the L subtree before revisiting this same
                            // child as I. Its descendants are searched again.
                            *tag.last_mut().expect("alternate has a residue") = b'I';
                            if tag.len() >= self.options.min_tag_length {
                                work.emit(&tag, tags)?;
                            }
                            frame.next = frame.index + 1;
                            frame.alternate = false;
                            stack.push(frame);
                        } else if !stack.is_empty() {
                            tag.pop();
                        }
                        continue;
                    }
                    stack.last_mut().expect("active DFS frame").next += 1;
                    if let Some(residue) = residue {
                        tag.push(residue);
                        if tag.len() >= self.options.min_tag_length {
                            work.emit(&tag, tags)?;
                        }
                        stack.push(Frame {
                            index: next,
                            next: next + 1,
                            alternate: residue == b'L',
                        });
                    }
                }
                debug_assert!(tag.is_empty());
            }
        }
        Ok(())
    }

    fn amino_acid(&self, mass: f64, work: &mut TaggerWork) -> Result<Option<u8>> {
        work.consume(1)?;
        if mass < self.min_gap || mass > self.max_gap {
            return Ok(None);
        }
        let tolerance = delta(self.options.tolerance, mass)?;
        let lower = finite(mass - tolerance, "mass tolerance lower bound")?;
        let mut begin = 0;
        let mut end = self.masses.len();
        while begin < end {
            work.consume(1)?;
            let middle = begin + (end - begin) / 2;
            if self.masses[middle].0 < lower {
                begin = middle + 1;
            } else {
                end = middle;
            }
        }
        if begin == self.masses.len() {
            return Ok(None);
        }
        let mut distance = finite((self.masses[begin].0 - mass).abs(), "residue mass error")?;
        // A mass exactly on the lower boundary rejects the whole query in the
        // source, even when a later table entry would be strictly inside.
        if distance >= tolerance {
            return Ok(None);
        }
        let mut best = begin;
        let mut best_distance = distance;
        while distance < tolerance {
            work.consume(1)?;
            begin += 1;
            if begin == self.masses.len() {
                break;
            }
            distance = finite((self.masses[begin].0 - mass).abs(), "residue mass error")?;
            if best_distance > distance {
                best = begin;
                best_distance = distance;
            }
        }
        Ok(Some(self.masses[best].1))
    }
}

#[derive(Clone, Copy)]
struct Frame {
    index: usize,
    next: usize,
    alternate: bool,
}

struct TaggerWork {
    remaining: usize,
    bytes: usize,
    text: usize,
    tags: usize,
}
impl Default for TaggerWork {
    fn default() -> Self {
        Self {
            remaining: MAX_TAGGER_WORK,
            bytes: MAX_TAGGER_BYTES,
            text: MAX_TAGGER_TEXT_BYTES,
            tags: 0,
        }
    }
}
impl TaggerWork {
    fn consume(&mut self, count: usize) -> Result<()> {
        self.remaining = self
            .remaining
            .checked_sub(count)
            .ok_or_else(|| invalid("tagger cumulative work limit exceeded"))?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize) -> Result<()> {
        let bytes = count
            .checked_mul(size_of::<T>())
            .ok_or_else(|| invalid("tagger allocation accounting overflow"))?;
        self.bytes = self
            .bytes
            .checked_sub(bytes)
            .ok_or_else(|| invalid("tagger cumulative allocation limit exceeded"))?;
        Ok(())
    }
    fn existing_tags(&mut self, tags: &[String]) -> Result<()> {
        if tags.len() > MAX_TAGGER_TAGS {
            return Err(invalid("tagger tag count limit exceeded"));
        }
        self.consume(tags.len())?;
        self.tags = tags.len();
        for tag in tags {
            self.text = self
                .text
                .checked_sub(tag.len())
                .ok_or_else(|| invalid("tagger cumulative text limit exceeded"))?;
        }
        Ok(())
    }
    fn emit(&mut self, tag: &[u8], tags: &mut Vec<String>) -> Result<()> {
        if self.tags == MAX_TAGGER_TAGS {
            return Err(invalid("tagger tag count limit exceeded"));
        }
        self.text = self
            .text
            .checked_sub(tag.len())
            .ok_or_else(|| invalid("tagger cumulative text limit exceeded"))?;
        self.consume(
            tag.len()
                .checked_mul(2)
                .ok_or_else(|| invalid("tagger text work overflow"))?,
        )?;
        self.allocate::<u8>(tag.len())?;
        if tags.len() == tags.capacity() {
            let capacity = tags.capacity().saturating_mul(2).clamp(1, MAX_TAGGER_TAGS);
            self.allocate::<String>(capacity)?;
            tags.reserve_exact(capacity - tags.len());
        }
        tags.push(
            String::from_utf8(tag.to_vec())
                .map_err(|_| invalid("invalid resolved tagger residue code"))?,
        );
        self.tags += 1;
        Ok(())
    }
    fn compare(&mut self, left: &str, right: &str) -> Result<Ordering> {
        for (left, right) in left.bytes().zip(right.bytes()) {
            self.consume(1)?;
            let order = left.cmp(&right);
            if order != Ordering::Equal {
                return Ok(order);
            }
        }
        self.consume(1)?;
        Ok(left.len().cmp(&right.len()))
    }
}

fn sorted_unique_indices(
    old: &[String],
    new: &[String],
    work: &mut TaggerWork,
) -> Result<Vec<usize>> {
    let count = old
        .len()
        .checked_add(new.len())
        .ok_or_else(|| invalid("tagger tag count overflow"))?;
    work.allocate::<usize>(count)?;
    work.allocate::<usize>(count)?;
    work.consume(count * 2)?;
    let at = |index: usize| {
        if index < old.len() {
            old[index].as_str()
        } else {
            new[index - old.len()].as_str()
        }
    };
    let mut order: Vec<_> = (0..count).collect();
    let mut buffer = vec![0; count];
    let mut width = 1;
    // Fallible bottom-up merge sort permits immediate budget failure. A std
    // comparator cannot return Result, and changing its order on failure could
    // violate sort's contract. Only bounded index vectors are staged here.
    while width < count {
        for start in (0..count).step_by(width * 2) {
            let middle = (start + width).min(count);
            let end = (middle + width).min(count);
            let (mut left, mut right) = (start, middle);
            for destination in &mut buffer[start..end] {
                work.consume(1)?;
                let take_left = right == end
                    || (left < middle
                        && work.compare(at(order[left]), at(order[right]))? != Ordering::Greater);
                *destination = if take_left {
                    let value = order[left];
                    left += 1;
                    value
                } else {
                    let value = order[right];
                    right += 1;
                    value
                };
            }
        }
        std::mem::swap(&mut order, &mut buffer);
        width *= 2;
    }
    let mut unique = 0;
    for read in 0..count {
        work.consume(1)?;
        if unique == 0 || work.compare(at(order[unique - 1]), at(order[read]))? != Ordering::Equal {
            order[unique] = order[read];
            unique += 1;
        }
    }
    order.truncate(unique);
    Ok(order)
}

fn delta(tolerance: Tolerance, mass: f64) -> Result<f64> {
    finite(
        match tolerance {
            Tolerance::Absolute(value) => value,
            Tolerance::Ppm(value) => (value / 1e6) * mass,
        },
        "mass tolerance",
    )
}
fn check_max_charge(charge: usize) -> Result<()> {
    if charge == usize::MAX {
        return Err(invalid(
            "tagger maximum charge would overflow the source inclusive loop",
        ));
    }
    Ok(())
}
fn check_peaks(count: usize) -> Result<()> {
    if count > MAX_TAGGER_PEAKS {
        return Err(invalid("tagger peak count limit exceeded"));
    }
    Ok(())
}
fn finite(value: f64, what: &str) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(format!("nonfinite tagger {what}")))
    }
}
fn invalid(message: impl Into<String>) -> Error {
    Error::InvalidValue(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deep_paths_use_heap_frames_and_branching_work_failure_is_atomic() {
        let mut options = TaggerOptions::new(9_999, Tolerance::Absolute(0.01));
        options.max_tag_length = 9_999;
        let mut tagger = Tagger {
            options,
            masses: vec![(1.0, b'A')],
            min_gap: 0.99,
            max_gap: 1.01,
        };
        let mzs: Vec<_> = (0..10_000).map(f64::from).collect();
        let tags = tagger.get_tags(&mzs).unwrap();
        assert_eq!(tags, ["A".repeat(9_999)]);

        tagger.options.min_tag_length = 1;
        tagger.options.max_tag_length = 20;
        tagger.masses[0].1 = b'L';
        let mut old = vec!["kept".to_owned()];
        let pointer = old[0].as_ptr();
        let mut work = TaggerWork {
            remaining: 500,
            ..Default::default()
        };
        assert!(
            tagger
                .append_with_work(&mzs[..31], &mut old, &mut work)
                .unwrap_err()
                .to_string()
                .contains("work limit")
        );
        assert!(
            work.tags > 1,
            "some DFS paths were emitted before exhaustion"
        );
        assert_eq!(old, ["kept"]);
        assert_eq!(old[0].as_ptr(), pointer);
    }

    #[test]
    fn sort_comparison_budget_failure_preserves_old_strings_and_their_allocations() {
        let tagger = Tagger::new(TaggerOptions::new(0, Tolerance::Absolute(0.01))).unwrap();
        let mut tags = vec!["a".repeat(1000) + "b", "a".repeat(1000) + "a"];
        let saved = tags.clone();
        let pointers: Vec<_> = tags.iter().map(|tag| tag.as_ptr()).collect();
        let mut work = TaggerWork {
            remaining: 100,
            ..Default::default()
        };
        assert!(
            tagger
                .append_with_work(&[], &mut tags, &mut work)
                .unwrap_err()
                .to_string()
                .contains("work limit")
        );
        assert_eq!(tags, saved);
        assert_eq!(
            tags.iter().map(|tag| tag.as_ptr()).collect::<Vec<_>>(),
            pointers
        );
    }

    #[test]
    fn source_lower_boundary_rejects_later_inside_mass_and_ties_keep_lower_mass() {
        let options = TaggerOptions::new(1, Tolerance::Absolute(2.0));
        let mut tagger = Tagger {
            options,
            masses: vec![(100.0, b'A'), (101.0, b'B'), (104.0, b'C')],
            min_gap: 98.0,
            max_gap: 106.0,
        };
        assert_eq!(
            tagger
                .amino_acid(102.0, &mut TaggerWork::default())
                .unwrap(),
            None
        );
        tagger.options.tolerance = Tolerance::Absolute(3.0);
        assert_eq!(
            tagger
                .amino_acid(102.5, &mut TaggerWork::default())
                .unwrap(),
            Some(b'B')
        );
    }
}
