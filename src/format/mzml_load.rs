// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! PeakFileOptions execution for the currently represented mzML data model.

use super::{Error, Record, Result};
use crate::format::peak_options::PeakFileOptions;
use crate::kernel::{DataArray, NumericRange};

/// Scientific loading choices alongside independent selection resource limits.
/// General XML/binary limits are still supplied as the existing ReadOptions.
#[derive(Clone, Debug)]
pub struct LoadOptions {
    pub scientific: PeakFileOptions,
    /// Native extension for chromatogram-only reads. Skipped records remain validated.
    pub skip_spectra: bool,
    /// Cumulative work for range/membership checks, sorting and aligned selection.
    pub max_selection_work: usize,
    /// Cumulative selection index storage; never includes or resets binary limits.
    pub max_selection_bytes: usize,
}
impl Default for LoadOptions {
    fn default() -> Self {
        Self {
            scientific: PeakFileOptions::default(),
            skip_spectra: false,
            max_selection_work: 500_000_000,
            max_selection_bytes: 256 * 1024 * 1024,
        }
    }
}
impl LoadOptions {
    pub(super) fn validate(&self) -> Result<()> {
        let unsupported = if self.scientific.skip_xml_checks {
            Some("disabling XML checks")
        } else if !self.scientific.precursor_mz_selected_ion {
            Some("isolation-target precursor selection")
        } else {
            None
        };
        if let Some(option) = unsupported {
            return Err(Error::Unsupported(format!(
                "mzML loading does not implement {option}"
            )));
        }
        Ok(())
    }
}
fn limit() -> Error {
    Error::InvalidValue("mzML selection resource limit exceeded".into())
}
fn product(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
/// Source DRange::encloses uses these two comparisons. This preserves its
/// half-open endpoints and NaN-bound behavior, unlike Range::contains.
pub(super) fn contains(range: NumericRange, value: f64) -> bool {
    !(value < range.min || value >= range.max)
}

pub(super) struct State<'a> {
    pub options: &'a LoadOptions,
    work: usize,
    bytes: usize,
}
impl<'a> State<'a> {
    pub fn new(options: &'a LoadOptions) -> Result<Self> {
        options.validate()?;
        let mut state = Self {
            options,
            work: options.max_selection_work,
            bytes: options.max_selection_bytes,
        };
        state.spend(1)?;
        Ok(state)
    }
    pub(super) fn spend(&mut self, amount: usize) -> Result<()> {
        self.work = self.work.checked_sub(amount).ok_or_else(limit)?;
        Ok(())
    }
    /// Selection happens while decoded f64 arrays are still available. The
    /// existing record validation then runs even for excluded whole records.
    pub fn apply(
        &mut self,
        record: &mut Record,
        positions: &mut Vec<f64>,
        intensities: &mut Vec<f64>,
    ) -> Result<bool> {
        self.spend(1)?;
        let settings = &self.options.scientific;
        let spectrum = record.spectrum.is_some();
        let keep_record = if let Some(s) = &record.spectrum {
            let mut keep = !self.options.skip_spectra && !record.precursor_outside;
            if record.seen_fields.contains("ms_level") && settings.has_ms_levels() {
                self.spend(settings.ms_levels().len())?;
                keep &=
                    i32::try_from(s.ms_level).is_ok_and(|level| settings.contains_ms_level(level));
            }
            if record.rt_seen && settings.has_rt_range() {
                keep &= contains(settings.rt_range(), s.rt);
            }
            keep
        } else {
            !settings.skip_chromatograms
        };
        if !keep_record {
            return Ok(false);
        }
        let n = positions.len();
        let (floats, integers, strings) = if let Some(s) = &mut record.spectrum {
            (
                &mut s.float_data_arrays,
                &mut s.integer_data_arrays,
                &mut s.string_data_arrays,
            )
        } else {
            let c = record.chromatogram.as_mut().unwrap();
            (
                &mut c.float_data_arrays,
                &mut c.integer_data_arrays,
                &mut c.string_data_arrays,
            )
        };
        let raw_floats = &mut record.raw_float_arrays;
        let arrays = floats
            .len()
            .checked_add(raw_floats.len())
            .and_then(|n| n.checked_add(integers.len()))
            .and_then(|n| n.checked_add(strings.len()))
            .ok_or_else(limit)?;
        self.spend(arrays)?;
        for length in floats
            .iter()
            .map(|a| a.data.len())
            .chain(raw_floats.iter().map(|a| a.data.len()))
            .chain(integers.iter().map(|a| a.data.len()))
            .chain(strings.iter().map(|a| a.data.len()))
        {
            if length != 0 && length != n {
                return Err(Error::InvalidValue(
                    "unaligned mzML auxiliary array before selection".into(),
                ));
            }
        }
        let coordinate_range = if spectrum {
            settings.has_mz_range().then(|| settings.mz_range())
        } else {
            settings.has_rt_range().then(|| settings.rt_range())
        };
        let intensity_range = settings
            .has_intensity_range()
            .then(|| settings.intensity_range());
        let sort = if spectrum {
            settings.sort_spectra_by_mz
        } else {
            settings.sort_chromatograms_by_rt
        };
        self.spend(n)?;
        let needs_sort = sort && positions.windows(2).any(|p| p[0] > p[1]);
        if coordinate_range.is_none() && intensity_range.is_none() && !needs_sort {
            return Ok(true);
        }
        // One selected-index vector plus the two inverse-permutation vectors.
        let scratch = product(product(n, 3)?, std::mem::size_of::<usize>())?;
        self.bytes = self.bytes.checked_sub(scratch).ok_or_else(limit)?;
        let visits = product(n, arrays.checked_add(8).ok_or_else(limit)?)?;
        self.spend(visits)?;
        let mut indices = Vec::with_capacity(n);
        for i in 0..n {
            if coordinate_range.is_none_or(|r| contains(r, positions[i]))
                && intensity_range.is_none_or(|r| contains(r, intensities[i]))
            {
                indices.push(i);
            }
        }
        if needs_sort {
            let levels = (usize::BITS - indices.len().max(1).leading_zeros()) as usize;
            // Conservative sort comparison allowance; no element/string cloning.
            self.spend(product(product(indices.len(), levels)?, 8)?)?;
            indices.sort_unstable_by(|&a, &b| {
                positions[a]
                    .partial_cmp(&positions[b])
                    .unwrap()
                    .then(a.cmp(&b))
            });
        }
        if !indices.is_empty() {
            let mut original_at: Vec<_> = (0..n).collect();
            let mut position_of: Vec<_> = (0..n).collect();
            for (to, &original) in indices.iter().enumerate() {
                let from = position_of[original];
                if from == to {
                    continue;
                }
                let displaced = original_at[to];
                positions.swap(to, from);
                intensities.swap(to, from);
                swap_arrays(floats, to, from);
                swap_arrays(raw_floats, to, from);
                swap_arrays(integers, to, from);
                swap_arrays(strings, to, from);
                original_at.swap(to, from);
                position_of[original] = to;
                position_of[displaced] = from;
            }
        }
        positions.truncate(indices.len());
        intensities.truncate(indices.len());
        truncate_arrays(floats, indices.len());
        truncate_arrays(raw_floats, indices.len());
        truncate_arrays(integers, indices.len());
        truncate_arrays(strings, indices.len());
        Ok(true)
    }
}
fn swap_arrays<T>(arrays: &mut [DataArray<T>], a: usize, b: usize) {
    for array in arrays {
        if !array.data.is_empty() {
            array.data.swap(a, b);
        }
    }
}
fn truncate_arrays<T>(arrays: &mut [DataArray<T>], length: usize) {
    for array in arrays {
        array.data.truncate(length);
    }
}
