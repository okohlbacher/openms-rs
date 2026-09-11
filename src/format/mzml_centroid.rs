// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! One-pass spectrum classification through the ordinary consumer.
use super::{LoadOptions, MSDataConsumer, ReadOptions, TransformOptions, transform_with_options};
use crate::{
    Error, MSChromatogram, MSSpectrum, Result, kernel::SpectrumType, metadata::ExperimentalSettings,
};
use std::{collections::BTreeMap, ops::ControlFlow, path::Path};

/// Classification/result-map budgets, shared across every delivered spectrum.
/// max_points applies per actual estimate; parser/consumer limits are separate.
pub use crate::kernel::SpectrumTypeQueryLimits as CentroidInfoLimits;
/// Source per-MS-level classification counters.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SpecInfo {
    pub count_centroided: usize,
    pub count_profile: usize,
    pub count_unknown: usize,
}
/// Inspect until ten recognized spectra have been counted across all MS levels.
pub fn centroid_info(path: impl AsRef<Path>) -> Result<BTreeMap<u32, SpecInfo>> {
    centroid_info_with_options(
        path,
        10,
        &LoadOptions::default(),
        &ReadOptions::default(),
        CentroidInfoLimits::default(),
    )
}
/// Unknown spectra are counted but do not consume the positive global quota.
/// Locally forces data population without changing caller options (CPP-018).
/// Zero quota is rejected before opening input (CPP-048). Pool predecoding and
/// scientific filters retain the consumer's normal behavior; no setup pass runs.
pub fn centroid_info_with_options(
    path: impl AsRef<Path>,
    first_n_spectra_only: usize,
    options: &LoadOptions,
    read: &ReadOptions,
    mut limits: CentroidInfoLimits,
) -> Result<BTreeMap<u32, SpecInfo>> {
    if first_n_spectra_only == 0 {
        return Err(Error::InvalidValue(
            "centroid inspection quota must be positive".into(),
        ));
    }
    // PeakFileOptions owns only its MS-level vector. Charge its clone before
    // creating local forced options; all other configuration fields are scalar.
    let n = options.scientific.ms_levels().len();
    let mut meter = crate::kernel::data_array::Meter {
        work: &mut limits.max_work,
        bytes: &mut limits.max_bytes,
    };
    meter.charge(1, 0)?;
    meter.slots::<i32>(n)?;
    let mut load = options.clone();
    load.scientific.fill_data = true;
    let transform = TransformOptions {
        skip_first_pass: true,
        skip_full_count: true,
        load,
        read: *read,
        ..Default::default()
    };
    let mut counter = Counter {
        remaining: first_n_spectra_only,
        limits,
        counts: BTreeMap::new(),
    };
    transform_with_options(path, &mut counter, &transform)?;
    Ok(counter.counts)
}
struct Counter {
    remaining: usize,
    limits: CentroidInfoLimits,
    counts: BTreeMap<u32, SpecInfo>,
}
impl MSDataConsumer for Counter {
    fn set_expected_size(&mut self, _: usize, _: usize) -> Result<()> {
        Err(Error::InvalidValue(
            "unexpected centroid setup callback".into(),
        ))
    }
    fn set_experimental_settings(&mut self, _: &ExperimentalSettings) -> Result<()> {
        Err(Error::InvalidValue(
            "unexpected centroid settings callback".into(),
        ))
    }
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        let kind = spectrum.get_type_with_budget(
            true,
            self.limits.max_points,
            &mut self.limits.max_work,
            &mut self.limits.max_bytes,
        )?;
        let level = spectrum.ms_level;
        let height = (usize::BITS - self.counts.len().leading_zeros()) as usize;
        let mut meter = crate::kernel::data_array::Meter {
            work: &mut self.limits.max_work,
            bytes: &mut self.limits.max_bytes,
        };
        // Two bounded-key BTree lookups (contains + entry), then one counter.
        meter.charge(2 + 24 * height, 0)?;
        if !self.counts.contains_key(&level) {
            meter.tree::<(u32, SpecInfo)>(1)?;
        }
        let entry = self.counts.entry(level).or_default();
        let count = match kind {
            SpectrumType::Centroid => &mut entry.count_centroided,
            SpectrumType::Profile => &mut entry.count_profile,
            SpectrumType::Unknown => &mut entry.count_unknown,
        };
        *count = count
            .checked_add(1)
            .ok_or_else(|| Error::InvalidValue("centroid count overflow".into()))?;
        if kind != SpectrumType::Unknown {
            self.remaining -= 1;
        }
        Ok(if self.remaining == 0 {
            ControlFlow::Break(())
        } else {
            ControlFlow::Continue(())
        })
    }
    fn consume_chromatogram(&mut self, _: &mut MSChromatogram) -> Result<ControlFlow<()>> {
        Ok(ControlFlow::Continue(()))
    }
}
