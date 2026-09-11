// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::MetaInfo;
use crate::{Error, Result};
use std::mem::size_of;

/// Complete source contact fields; ordinary owned fields replace accessors.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContactPerson {
    pub first_name: String,
    pub last_name: String,
    pub institution: String,
    pub email: String,
    pub contact_info: String,
    pub url: String,
    pub address: String,
    pub metadata: MetaInfo,
}
impl ContactPerson {
    /// Source name splitting: comma wins and trims its first two fields; the
    /// space branch retains empty fields. Further fields are ignored. A name
    /// without either delimiter changes only last_name, retaining first_name.
    pub fn set_name(&mut self, name: &str) -> Result<()> {
        if name.len() > 4 * 1024 * 1024 {
            return Err(bad("contact name exceeds 4 MiB"));
        }
        let pair = if name.contains(',') {
            let mut parts = name.split(',');
            let last = crate::data_structures::list::trim(parts.next().unwrap());
            let first = crate::data_structures::list::trim(parts.next().unwrap());
            Some((first, last))
        } else if name.contains(' ') {
            let mut parts = name.split(' ');
            Some((parts.next().unwrap(), parts.next().unwrap()))
        } else {
            None
        };
        if let Some((first, last)) = pair {
            let first = copy_string(first)?;
            let last = copy_string(last)?;
            self.first_name = first;
            self.last_name = last;
        } else {
            self.last_name = copy_string(name)?;
        }
        Ok(())
    }
}

/// Source chromatography instrument/settings and an owned gradient.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HPLC {
    pub instrument: String,
    pub column: String,
    /// Degrees Celsius, source default21.
    pub temperature: i32,
    /// Bar; native unsigned storage matches the public source API.
    pub pressure: u32,
    /// Microliters per second.
    pub flux: u32,
    pub comment: String,
    pub gradient: Gradient,
}
impl Default for HPLC {
    fn default() -> Self {
        Self {
            instrument: String::new(),
            column: String::new(),
            temperature: 21,
            pressure: 0,
            flux: 0,
            comment: String::new(),
            gradient: Gradient::default(),
        }
    }
}

/// Named eluents, strictly increasing signed timepoints, and row-major integer
/// percentages. Source clear-axis operations deliberately retain old percentage
/// storage; clear_percentages explicitly rebuilds its shape with zeros.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Gradient {
    eluents: Vec<String>,
    times: Vec<i32>,
    percentages: Vec<Vec<u32>>,
}
impl Gradient {
    /// Aggregate axes, row descriptors and stored percentage cells.
    pub const MAX_ITEMS: usize = 1_000_000;
    /// Conservative logical storage including growth/scratch headroom.
    pub const MAX_BYTES: usize = 64 * 1024 * 1024;

    pub fn new() -> Self {
        Self::default()
    }
    pub fn eluents(&self) -> &[String] {
        &self.eluents
    }
    pub fn timepoints(&self) -> &[i32] {
        &self.times
    }
    pub fn percentages(&self) -> &[Vec<u32>] {
        &self.percentages
    }
    pub fn add_eluent(&mut self, eluent: &str) -> Result<()> {
        if self.eluents.iter().any(|item| item == eluent) {
            return Err(bad("gradient eluent already exists"));
        }
        let extra = add(2, self.times.len())?;
        let bytes = add(
            add(size_of::<String>(), size_of::<Vec<u32>>())?,
            add(eluent.len(), mul(self.times.len(), size_of::<u32>())?)?,
        )?;
        self.preflight(extra, bytes)?;
        let name = copy_string(eluent)?;
        let mut row = Vec::new();
        reserve(&mut row, self.times.len())?;
        row.resize(self.times.len(), 0);
        reserve(&mut self.eluents, 1)?;
        reserve(&mut self.percentages, 1)?;
        self.eluents.push(name);
        self.percentages.push(row);
        Ok(())
    }
    pub fn add_timepoint(&mut self, timepoint: i32) -> Result<()> {
        if self.times.last().is_some_and(|last| timepoint <= *last) {
            return Err(bad("gradient timepoints must strictly increase"));
        }
        self.preflight(
            add(1, self.eluents.len())?,
            mul(add(1, self.eluents.len())?, 4)?,
        )?;
        if self.percentages.len() < self.eluents.len() {
            return Err(bad("gradient percentage storage has too few rows"));
        }
        reserve(&mut self.times, 1)?;
        for row in &mut self.percentages[..self.eluents.len()] {
            reserve(row, 1)?;
        }
        self.times.push(timepoint);
        for row in &mut self.percentages[..self.eluents.len()] {
            row.push(0);
        }
        Ok(())
    }
    /// Clears only the names, preserving the source's stale percentage rows.
    pub fn clear_eluents(&mut self) {
        self.eluents.clear();
    }
    /// Clears only the axis, preserving the source's stale percentage columns.
    pub fn clear_timepoints(&mut self) {
        self.times.clear();
    }
    pub fn percentage(&self, eluent: &str, timepoint: i32) -> Result<u32> {
        let (row, column) = self.indices(eluent, timepoint)?;
        self.percentages
            .get(row)
            .and_then(|row| row.get(column))
            .copied()
            .ok_or_else(|| bad("gradient percentage storage does not cover its axes"))
    }
    pub fn set_percentage(&mut self, eluent: &str, timepoint: i32, percentage: u32) -> Result<()> {
        let (row, column) = self.indices(eluent, timepoint)?;
        if percentage > 100 {
            return Err(bad("gradient percentage exceeds 100"));
        }
        let value = self
            .percentages
            .get_mut(row)
            .and_then(|row| row.get_mut(column))
            .ok_or_else(|| bad("gradient percentage storage does not cover its axes"))?;
        *value = percentage;
        Ok(())
    }
    /// Rebuilds rows to current axes with all percentages zero, repairing any
    /// stale shape left by the source clear-axis operations. Atomic on error.
    pub fn clear_percentages(&mut self) -> Result<()> {
        let cells = mul(self.eluents.len(), self.times.len())?;
        let descriptors = add(self.eluents.len(), cells)?;
        let bytes = add(
            mul(self.eluents.len(), size_of::<Vec<u32>>())?,
            mul(cells, 4)?,
        )?;
        // Include the old table, which remains alive until replacement succeeds.
        self.preflight(descriptors, bytes)?;
        let mut rows = Vec::new();
        reserve(&mut rows, self.eluents.len())?;
        for _ in &self.eluents {
            let mut row = Vec::new();
            reserve(&mut row, self.times.len())?;
            row.resize(self.times.len(), 0);
            rows.push(row);
        }
        self.percentages = rows;
        Ok(())
    }
    /// All active timepoint columns must total100. No active timepoints is
    /// valid; active times without eluents is invalid. Undefined source indexing
    /// from stale shape becomes a checked error.
    pub fn is_valid(&self) -> Result<bool> {
        for column in 0..self.times.len() {
            let mut total = 0u64;
            for row in 0..self.eluents.len() {
                let value = self
                    .percentages
                    .get(row)
                    .and_then(|row| row.get(column))
                    .ok_or_else(|| bad("gradient percentage storage does not cover its axes"))?;
                total += u64::from(*value);
            }
            if total != 100 {
                return Ok(false);
            }
        }
        Ok(true)
    }
    fn indices(&self, eluent: &str, timepoint: i32) -> Result<(usize, usize)> {
        let row = self
            .eluents
            .iter()
            .position(|value| value == eluent)
            .ok_or_else(|| bad("unknown gradient eluent"))?;
        let column = self
            .times
            .binary_search(&timepoint)
            .map_err(|_| bad("unknown gradient timepoint"))?;
        Ok((row, column))
    }
    fn preflight(&self, extra_items: usize, extra_bytes: usize) -> Result<()> {
        let mut items = add(
            add(self.eluents.len(), self.times.len())?,
            self.percentages.len(),
        )?;
        let mut bytes = add(
            mul(self.eluents.capacity(), size_of::<String>())?,
            mul(self.times.capacity(), 4)?,
        )?;
        bytes = add(
            bytes,
            mul(self.percentages.capacity(), size_of::<Vec<u32>>())?,
        )?;
        for name in &self.eluents {
            bytes = add(bytes, name.capacity())?;
        }
        for row in &self.percentages {
            items = add(items, row.len())?;
            bytes = add(bytes, mul(row.capacity(), 4)?)?;
        }
        if add(items, extra_items)? > Self::MAX_ITEMS
            || mul(add(bytes, extra_bytes)?, 4)? > Self::MAX_BYTES
        {
            return Err(bad("gradient exceeds descriptor or byte limit"));
        }
        Ok(())
    }
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.into())
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b)
        .ok_or_else(|| bad("gradient size overflow"))
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b)
        .ok_or_else(|| bad("gradient size overflow"))
}
fn reserve<T>(values: &mut Vec<T>, additional: usize) -> Result<()> {
    let wanted = add(values.len(), additional)?;
    if wanted > values.capacity() {
        // Retain source amortized append cost. The fourfold preflight storage
        // allowance covers geometric capacity and simultaneous old/new buffers.
        let capacity = wanted.max(values.capacity().saturating_mul(2)).max(4);
        values
            .try_reserve_exact(capacity - values.len())
            .map_err(|_| bad("cannot allocate gradient values"))?;
    }
    Ok(())
}
fn copy_string(text: &str) -> Result<String> {
    let mut value = String::new();
    value
        .try_reserve_exact(text.len())
        .map_err(|_| bad("cannot allocate metadata text"))?;
    value.push_str(text);
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn empty_axes_retained_capacity_is_charged_before_growth() {
        let capacity = Gradient::MAX_BYTES / (4 * size_of::<String>());
        let mut gradient = Gradient {
            eluents: Vec::with_capacity(capacity),
            ..Default::default()
        };
        assert!(gradient.add_eluent("new").is_err());
        assert!(gradient.eluents.is_empty() && gradient.percentages.is_empty());
        assert_eq!(gradient.eluents.capacity(), capacity);
    }
}
