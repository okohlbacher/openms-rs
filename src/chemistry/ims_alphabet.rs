// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Owned IMS alphabet and replaceable plain-text parser, from OpenMS 54a232f.

use super::{IMSElement, IMSIsotopeOptions, IMSIsotopePeak};
use crate::{Error, Result};
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

#[path = "ims_alphabet_parser.rs"]
mod parser;
pub use parser::{IMSAlphabetParser, IMSAlphabetTextParser};

const MAX_ELEMENTS: usize = 100_000;
const MAX_LABEL: usize = 1024 * 1024;
const MAX_BYTES: usize = 64 * 1024 * 1024;
const MAX_TEXT: usize = 8 * 1024 * 1024;
const MAX_WORK: usize = 50_000_000;

/// An indexed alphabet, retaining insertion order and duplicate names.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct IMSAlphabet {
    elements: Vec<IMSElement>,
    payload_bytes: usize,
}

impl IMSAlphabet {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_elements(elements: Vec<IMSElement>) -> Result<Self> {
        let mut work = Work::default();
        check_count(elements.len())?;
        work.consume(elements.len())?;
        let mut payload_bytes = 0usize;
        for element in &elements {
            payload_bytes = add(payload_bytes, element_bytes(element)?)?;
            check_payload(payload_bytes)?;
        }
        Ok(Self {
            elements,
            payload_bytes,
        })
    }

    pub fn len(&self) -> usize {
        self.elements.len()
    }
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty()
    }
    pub fn elements(&self) -> &[IMSElement] {
        &self.elements
    }

    pub fn element(&self, index: usize) -> Result<&IMSElement> {
        self.elements
            .get(index)
            .ok_or_else(|| invalid("IMS alphabet index is out of range"))
    }

    /// First matching name, using the source linear scan.
    pub fn get(&self, name: &str) -> Result<&IMSElement> {
        let index = self
            .find(name, &mut Work::default())?
            .ok_or_else(|| invalid("name was not found in IMS alphabet"))?;
        Ok(&self.elements[index])
    }

    pub fn name(&self, index: usize) -> Result<&str> {
        Ok(self.element(index)?.name())
    }
    pub fn mass(&self, index: usize) -> Result<f64> {
        self.element(index)?.mass(0)
    }
    pub fn mass_by_name(&self, name: &str) -> Result<f64> {
        self.get(name)?.mass(0)
    }
    pub fn has_name(&self, name: &str) -> Result<bool> {
        Ok(self.find(name, &mut Work::default())?.is_some())
    }

    pub fn masses(&self, isotope_index: usize) -> Result<Vec<f64>> {
        let mut result = Work::default().vector(self.len())?;
        for element in &self.elements {
            result.push(element.mass(isotope_index)?);
        }
        Ok(result)
    }

    pub fn average_masses(&self) -> Result<Vec<f64>> {
        let mut work = Work::default();
        let mut result = work.vector(self.len())?;
        for element in &self.elements {
            work.consume(mul(element.isotope_distribution().stored_len(), 4)?)?;
            result.push(element.average_mass()?);
        }
        Ok(result)
    }

    /// Replaces only the first match with a fresh name/mass element. If absent,
    /// `forced=false` is a no-op, even when the unused mass is nonfinite.
    pub fn set_element(&mut self, name: &str, mass: f64, forced: bool) -> Result<()> {
        let mut work = Work::default();
        let found = self.find(name, &mut work)?;
        if found.is_none() && !forced {
            return Ok(());
        }
        work.copy(mul(name.len(), 2)?)?;
        let replacement = IMSElement::from_mass(name, mass)?;
        if let Some(index) = found {
            let payload = add(
                self.payload_bytes - element_bytes(&self.elements[index])?,
                element_bytes(&replacement)?,
            )?;
            check_payload(payload)?;
            self.elements[index] = replacement;
            self.payload_bytes = payload;
        } else {
            self.push_with_work(replacement, &mut work)?;
        }
        Ok(())
    }

    pub fn push(&mut self, element: IMSElement) -> Result<()> {
        self.push_with_work(element, &mut Work::default())
    }

    pub fn push_mass(&mut self, name: &str, mass: f64) -> Result<()> {
        check_label(name)?;
        let mut work = Work::default();
        work.copy(mul(name.len(), 2)?)?;
        self.push_with_work(IMSElement::from_mass(name, mass)?, &mut work)
    }

    /// Removes only the first matching element.
    pub fn erase(&mut self, name: &str) -> Result<bool> {
        let mut work = Work::default();
        let Some(index) = self.find(name, &mut work)? else {
            return Ok(false);
        };
        work.consume(self.len() - index)?;
        let bytes = element_bytes(&self.elements[index])?;
        self.elements.remove(index);
        self.payload_bytes -= bytes;
        Ok(true)
    }

    /// Normal Rust destruction, retaining the allocation for the element vector.
    pub fn clear(&mut self) {
        self.elements.clear();
        self.payload_bytes = 0;
    }

    /// Stable ties are a deterministic native replacement for unspecified C++ ties.
    pub fn sort_by_names(&mut self) -> Result<()> {
        let mut work = Work::default();
        let order = sorted_indices(self.len(), &mut work, |a, b, work| {
            let a = self.elements[a].name();
            let b = self.elements[b].name();
            work.consume(add(a.len().min(b.len()), 1)?)?;
            Ok(a.cmp(b))
        })?;
        self.apply_order(order, &mut work)
    }

    pub fn sort_by_mass(&mut self) -> Result<()> {
        self.sort_mass_with_work(&mut Work::default())
    }

    /// Default plain-text parsing and mass sorting share one native budget.
    pub fn read(reader: impl BufRead) -> Result<Self> {
        let mut work = Work::default();
        let parsed = parser::parse_map(reader, &mut work)?;
        Self::from_map(&parsed, &mut work)
    }

    /// Source plain-file replacement; no compression or invented source metadata.
    pub fn load(&mut self, path: impl AsRef<Path>) -> Result<()> {
        let reader = BufReader::new(std::fs::File::open(path)?);
        let mut work = Work::default();
        work.consume(mul(self.len(), 4)?)?; // old element/string/vector destruction
        let parsed = parser::parse_map(reader, &mut work)?;
        *self = Self::from_map(&parsed, &mut work)?;
        Ok(())
    }

    /// The alphabet is atomic; the supplied parser owns its own processing budget
    /// and side effects. Its entire returned map is checked before alphabet commit.
    pub fn load_with_parser(
        &mut self,
        path: impl AsRef<Path>,
        parser: &mut dyn IMSAlphabetParser,
    ) -> Result<()> {
        let mut work = Work::default();
        work.consume(mul(self.len(), 4)?)?;
        parser.load(path.as_ref())?;
        *self = Self::from_map(parser.elements(), &mut work)?;
        Ok(())
    }

    /// Source verbose stream output. It is not a codec for the flat input format.
    pub fn to_text(&self, options: IMSIsotopeOptions) -> Result<String> {
        let mut work = Work::default();
        work.consume(self.len())?;
        let mut capacity = 0usize;
        for element in &self.elements {
            let bins = element.isotope_distribution().size(options)?;
            let labels = add(
                add(element.name().len(), element.sequence().len())?,
                "name:\t\nsequence:\t\nisotope distribution:\n\n".len(),
            )?;
            let text = add(add(labels, mul(bins, 32)?)?, 1)?;
            capacity = add(capacity, text)?;
            if capacity > MAX_TEXT {
                return Err(invalid("IMS alphabet text exceeds 8 MiB"));
            }
            // Account the frozen element/distribution formatter's scalar scratch,
            // intermediate strings and this final destination before invoking it.
            work.copy(add(mul(bins, 1024)?, mul(text, 3)?)?)?;
        }
        let mut result = String::new();
        result.try_reserve_exact(capacity).map_err(|_| limit())?;
        for element in &self.elements {
            result.push_str(&element.to_text(options)?);
            result.push('\n');
        }
        Ok(result)
    }

    /// Validates and formats the whole alphabet before writing any bytes.
    pub fn write(&self, mut writer: impl Write, options: IMSIsotopeOptions) -> Result<()> {
        writer.write_all(self.to_text(options)?.as_bytes())?;
        Ok(())
    }

    fn find(&self, name: &str, work: &mut Work) -> Result<Option<usize>> {
        check_label(name)?;
        for (index, element) in self.elements.iter().enumerate() {
            work.consume(add(element.name().len().min(name.len()), 1)?)?;
            if element.name() == name {
                return Ok(Some(index));
            }
        }
        Ok(None)
    }

    fn push_with_work(&mut self, element: IMSElement, work: &mut Work) -> Result<()> {
        check_count(add(self.len(), 1)?)?;
        let payload = add(self.payload_bytes, element_bytes(&element)?)?;
        check_payload(payload)?;
        if self.elements.len() == self.elements.capacity() {
            work.copy(mul(add(self.len(), 1)?, std::mem::size_of::<IMSElement>())?)?;
            self.elements.try_reserve_exact(1).map_err(|_| limit())?;
        }
        self.elements.push(element);
        self.payload_bytes = payload;
        Ok(())
    }

    fn from_map(elements: &BTreeMap<String, f64>, work: &mut Work) -> Result<Self> {
        parser::measure_map(elements, work)?;
        let mut output = Self {
            elements: work.vector(elements.len())?,
            payload_bytes: 0,
        };
        for (name, &mass) in elements {
            work.copy(add(
                mul(name.len(), 2)?,
                std::mem::size_of::<IMSIsotopePeak>(),
            )?)?;
            let element = IMSElement::from_mass(name, mass)?;
            output.payload_bytes = add(output.payload_bytes, element_bytes(&element)?)?;
            check_payload(output.payload_bytes)?;
            output.elements.push(element);
        }
        output.sort_mass_with_work(work)?;
        Ok(output)
    }

    fn sort_mass_with_work(&mut self, work: &mut Work) -> Result<()> {
        // Source std::sort never requests a mass for zero or one element.
        if self.len() < 2 {
            return Ok(());
        }
        let mut masses = work.vector(self.len())?;
        for element in &self.elements {
            masses.push(element.mass(0)?);
        }
        let order = sorted_indices(self.len(), work, |a, b, work| {
            work.consume(1)?;
            // All masses are finite; ordinary comparison equates signed zero.
            Ok(masses[a].partial_cmp(&masses[b]).unwrap())
        })?;
        self.apply_order(order, work)
    }

    fn apply_order(&mut self, order: Vec<usize>, work: &mut Work) -> Result<()> {
        let mut destination = work.vector(self.len())?;
        destination.resize(self.len(), 0);
        work.consume(mul(self.len(), 3)?)?;
        for (new, &old) in order.iter().enumerate() {
            destination[old] = new;
        }
        // Every fallible step precedes this move-only permutation.
        for index in 0..self.len() {
            while destination[index] != index {
                let other = destination[index];
                self.elements.swap(index, other);
                destination.swap(index, other);
            }
        }
        Ok(())
    }
}

fn sorted_indices(
    len: usize,
    work: &mut Work,
    mut compare: impl FnMut(usize, usize, &mut Work) -> Result<Ordering>,
) -> Result<Vec<usize>> {
    let mut order = work.vector(len)?;
    order.extend(0..len);
    let mut scratch = work.vector(len)?;
    scratch.resize(len, 0);
    let mut width = 1;
    while width < len {
        let mut start = 0;
        while start < len {
            let middle = (start + width).min(len);
            let end = (middle + width).min(len);
            let (mut left, mut right) = (start, middle);
            for slot in &mut scratch[start..end] {
                work.consume(1)?;
                if right == end
                    || (left < middle
                        && compare(order[left], order[right], work)? != Ordering::Greater)
                {
                    *slot = order[left];
                    left += 1;
                } else {
                    *slot = order[right];
                    right += 1;
                }
            }
            start = end;
        }
        std::mem::swap(&mut order, &mut scratch);
        width *= 2;
    }
    Ok(order)
}

fn element_bytes(element: &IMSElement) -> Result<usize> {
    add(
        add(
            std::mem::size_of::<IMSElement>(),
            add(element.name().len(), element.sequence().len())?,
        )?,
        mul(
            element.isotope_distribution().stored_len(),
            std::mem::size_of::<IMSIsotopePeak>(),
        )?,
    )
}
fn check_count(count: usize) -> Result<()> {
    if count > MAX_ELEMENTS {
        Err(invalid("IMS alphabet exceeds 100,000 elements"))
    } else {
        Ok(())
    }
}
fn check_label(name: &str) -> Result<()> {
    if name.len() > MAX_LABEL {
        Err(invalid("IMS alphabet name exceeds 1 MiB"))
    } else {
        Ok(())
    }
}
fn check_payload(bytes: usize) -> Result<()> {
    if bytes > MAX_BYTES {
        Err(invalid("IMS alphabet payload exceeds 64 MiB"))
    } else {
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
fn invalid(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn limit() -> Error {
    invalid("IMS alphabet resource limit exceeded")
}

struct Work {
    remaining: usize,
    bytes: usize,
}
impl Default for Work {
    fn default() -> Self {
        Self {
            remaining: MAX_WORK,
            bytes: MAX_BYTES,
        }
    }
}
impl Work {
    fn consume(&mut self, amount: usize) -> Result<()> {
        self.remaining = self.remaining.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    fn allocate(&mut self, amount: usize) -> Result<()> {
        self.bytes = self.bytes.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    fn copy(&mut self, amount: usize) -> Result<()> {
        self.consume(amount)?;
        self.allocate(amount)
    }
    fn vector<T>(&mut self, len: usize) -> Result<Vec<T>> {
        self.consume(len)?;
        self.allocate(mul(len, std::mem::size_of::<T>())?)?;
        let mut values = Vec::new();
        values.try_reserve_exact(len).map_err(|_| limit())?;
        Ok(values)
    }
}
