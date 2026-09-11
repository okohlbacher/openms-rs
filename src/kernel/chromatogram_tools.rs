// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

use super::{
    ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor,
    acquisition_fields as fields, data_array::Meter,
};
use crate::metadata::{ChromatogramType, Product, ScanMode};
use crate::{Error, Result};
use std::{
    cmp::Ordering,
    collections::{BTreeMap, btree_map::Entry},
    mem::size_of,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Per-call ceilings for a spectrum/chromatogram conversion.
///
/// Native bounds with no source counterpart, checked before any output is
/// produced.
pub struct ChromatogramConversionLimits {
    pub max_input_records: usize,
    /// Final destination length, including its retained prefix.
    pub max_output_records: usize,
    pub max_points: usize,
    pub max_work: usize,
    pub max_bytes: usize,
}
impl Default for ChromatogramConversionLimits {
    fn default() -> Self {
        Self {
            max_input_records: 1_000_000,
            max_output_records: 1_000_000,
            max_points: 10_000_000,
            max_work: 50_000_000,
            max_bytes: 256 * 1024 * 1024,
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
/// Conversion between chromatograms and spectra.
///
/// Ports `OpenMS/KERNEL/ChromatogramTools.h`. The source class is a stateless
/// collection of template functions; this is a unit struct with associated
/// functions for the same reason.
pub struct ChromatogramTools {
    pub limits: ChromatogramConversionLimits,
}
/// Removed records are returned unchanged, avoiding an unbounded destructor of
/// unrelated annotations inside the checked conversion. Ordinary caller drop applies.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ChromatogramConversionReport {
    pub added_spectra: usize,
    pub added_chromatograms: usize,
    /// Selected SRM records that lacked one precursor and at least one peak.
    pub skipped_spectra: usize,
    pub removed_spectra: Vec<MSSpectrum>,
    pub removed_chromatograms: Vec<MSChromatogram>,
}
impl ChromatogramTools {
    /// Append one MS2 scan per chromatogram point, then remove all chromatograms.
    /// Copied acquisition/precursor/product payload is metered before cloning.
    pub fn convert_chromatograms_to_spectra(
        &self,
        experiment: &mut MSExperiment,
    ) -> Result<ChromatogramConversionReport> {
        let mut w = Work::new(self.limits);
        cap(
            experiment.chromatograms.len(),
            self.limits.max_input_records,
        )?;
        w.spend(experiment.chromatograms.len())?;
        let mut count = 0;
        for c in &experiment.chromatograms {
            count = add(count, c.len())?;
            cap(count, self.limits.max_points)?;
        }
        if count == 0 {
            return Ok(ChromatogramConversionReport {
                removed_chromatograms: std::mem::take(&mut experiment.chromatograms),
                ..Default::default()
            });
        }
        // Output preparation visits every input descriptor again, including
        // empty chromatograms. The zero-point shortcut above does not.
        w.spend(experiment.chromatograms.len())?;
        let total = add(experiment.spectra.len(), count)?;
        cap(total, self.limits.max_output_records)?;
        w.spend(total)?;
        let mut output = w.vector::<MSSpectrum>(total)?;
        let mut added = w.vector::<MSSpectrum>(count)?;
        for c in &experiment.chromatograms {
            for p in &c.peaks {
                w.spend(3)?;
                finite(p.rt)?;
                finite(f64::from(p.intensity))?;
                finite(c.product.mz)?;
                let mut meter = w.meter();
                fields::instrument(&mut meter, &c.instrument_settings)?;
                fields::acquisition(&mut meter, &c.acquisition_info)?;
                fields::source(&mut meter, &c.source_file)?;
                fields::precursor(&mut meter, &c.precursor)?;
                fields::product(&mut meter, &c.product)?;
                let mut peaks = w.vector(1)?;
                peaks.push(Peak1D::new(c.product.mz, p.intensity));
                let mut precursors = w.vector(1)?;
                precursors.push(c.precursor.clone());
                let mut products = w.vector(1)?;
                products.push(c.product.clone());
                let mut s = MSSpectrum {
                    peaks,
                    rt: p.rt,
                    ms_level: 2,
                    instrument_settings: c.instrument_settings.clone(),
                    acquisition_info: c.acquisition_info.clone(),
                    source_file: c.source_file.clone(),
                    precursors,
                    products,
                    ..Default::default()
                };
                match c.chromatogram_type {
                    ChromatogramType::SelectedReactionMonitoring => {
                        s.instrument_settings.scan_mode = ScanMode::SelectedReactionMonitoring
                    }
                    ChromatogramType::SelectedIonMonitoring => {
                        s.instrument_settings.scan_mode = ScanMode::SelectedIonMonitoring
                    }
                    _ => {}
                }
                added.push(s);
            }
        }
        // No fallible operation follows ownership transfer.
        output.append(&mut experiment.spectra);
        output.append(&mut added);
        experiment.spectra = output;
        Ok(ChromatogramConversionReport {
            added_spectra: count,
            removed_chromatograms: std::mem::take(&mut experiment.chromatograms),
            ..Default::default()
        })
    }
    /// Append XIC then SRM groups in ascending exact m/z key order. Points keep
    /// encounter order. Removal follows source SRM scan mode, even for skipped
    /// records; forced non-SRM records remain in the experiment.
    pub fn convert_spectra_to_chromatograms(
        &self,
        experiment: &mut MSExperiment,
        remove_spectra: bool,
        force_conversion: bool,
    ) -> Result<ChromatogramConversionReport> {
        let mut w = Work::new(self.limits);
        cap(experiment.spectra.len(), self.limits.max_input_records)?;
        w.spend(experiment.spectra.len())?;
        let mut groups = BTreeMap::<GroupKey, Group>::new();
        let mut points = 0usize;
        let mut skipped = 0usize;
        let mut removed = 0usize;
        for (index, s) in experiment.spectra.iter().enumerate() {
            let srm = s.instrument_settings.scan_mode == ScanMode::SelectedReactionMonitoring;
            if remove_spectra && srm {
                removed = add(removed, 1)?;
            }
            if !srm && !force_conversion {
                continue;
            }
            let one = s.precursors.len() == 1 && !s.peaks.is_empty();
            if !one && !force_conversion {
                skipped = add(skipped, 1)?;
                continue;
            }
            points = add(points, s.len())?;
            cap(points, self.limits.max_points)?;
            w.spend(mul(s.len(), 3)?)?;
            if !s.peaks.is_empty() {
                finite(s.rt)?;
            }
            let precursor = if one {
                finite(s.precursors[0].mz)?;
                Some(Key(s.precursors[0].mz))
            } else {
                None
            };
            for p in &s.peaks {
                finite(p.mz)?;
                finite(f64::from(p.intensity))?;
                let key = if let Some(q) = precursor {
                    GroupKey::Srm(q, Key(p.mz))
                } else {
                    GroupKey::Xic(Key(p.mz))
                };
                // Bound tree comparisons before lookup. Finite scalar keys avoid
                // the undefined ordering of source NaNs; +/-zero compare equal.
                w.spend(mul(
                    16,
                    add(groups.len().saturating_add(1).ilog2() as usize, 1)?,
                )?)?;
                let next_count = add(groups.len(), 1)?;
                let group = match groups.entry(key) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => {
                        cap(
                            add(experiment.chromatograms.len(), next_count)?,
                            self.limits.max_output_records,
                        )?;
                        w.meter().tree::<(GroupKey, Group)>(1)?;
                        entry.insert(Group {
                            first: index,
                            points: Vec::new(),
                        })
                    }
                };
                w.push(&mut group.points, ChromatogramPeak::new(s.rt, p.intensity))?;
            }
        }
        let added_count = groups.len();
        if added_count == 0 && removed == 0 {
            return Ok(ChromatogramConversionReport {
                skipped_spectra: skipped,
                ..Default::default()
            });
        }
        let total = add(experiment.chromatograms.len(), added_count)?;
        if added_count != 0 {
            cap(total, self.limits.max_output_records)?;
        }
        w.spend(add(total, experiment.spectra.len())?)?;
        let mut output = if added_count != 0 {
            w.vector::<MSChromatogram>(total)?
        } else {
            Vec::new()
        };
        let mut added = w.vector::<MSChromatogram>(added_count)?;
        let mut retained = if remove_spectra {
            w.vector::<MSSpectrum>(experiment.spectra.len() - removed)?
        } else {
            Vec::new()
        };
        let mut removed_spectra = w.vector::<MSSpectrum>(removed)?;
        for (key, group) in groups {
            let s = &experiment.spectra[group.first];
            let mut meter = w.meter();
            fields::instrument(&mut meter, &s.instrument_settings)?;
            fields::acquisition(&mut meter, &s.acquisition_info)?;
            fields::source(&mut meter, &s.source_file)?;
            let (precursor, product, chromatogram_type, native_id) = match key {
                GroupKey::Srm(_, product) => {
                    fields::precursor(&mut meter, &s.precursors[0])?;
                    meter.text(&s.native_id)?;
                    meter.charge(13, 13)?;
                    (
                        s.precursors[0].clone(),
                        Product {
                            mz: product.0,
                            ..Default::default()
                        },
                        ChromatogramType::SelectedReactionMonitoring,
                        format!("chromatogram={}", s.native_id),
                    )
                }
                GroupKey::Xic(mz) => {
                    meter.charge(256, 256)?;
                    meter.tree::<(String, crate::metadata::MetaValue)>(1)?;
                    let mut precursor = Precursor::new(mz.0, 0);
                    precursor.cv_terms.metadata.insert(
                        "description".into(),
                        format!("XIC @ {}", crate::param::value::format_float(mz.0, true)).into(),
                    );
                    (
                        precursor,
                        Product::default(),
                        ChromatogramType::Mass,
                        String::new(),
                    )
                }
            };
            added.push(MSChromatogram {
                peaks: group.points,
                precursor,
                product,
                chromatogram_type,
                native_id,
                instrument_settings: s.instrument_settings.clone(),
                acquisition_info: s.acquisition_info.clone(),
                source_file: s.source_file.clone(),
                ..Default::default()
            });
        }
        if added_count != 0 {
            output.append(&mut experiment.chromatograms);
            output.append(&mut added);
            experiment.chromatograms = output;
        }
        if remove_spectra {
            for s in std::mem::take(&mut experiment.spectra) {
                if s.instrument_settings.scan_mode == ScanMode::SelectedReactionMonitoring {
                    removed_spectra.push(s);
                } else {
                    retained.push(s);
                }
            }
            experiment.spectra = retained;
        }
        Ok(ChromatogramConversionReport {
            added_chromatograms: added_count,
            skipped_spectra: skipped,
            removed_spectra,
            ..Default::default()
        })
    }
}
// Keys are constructed only after finite validation. Exact numeric ordering
// intentionally equates signed zero, as std::map<double> does.
#[derive(Clone, Copy, Debug)]
struct Key(f64);
impl PartialEq for Key {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl Eq for Key {}
impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Key {
    fn cmp(&self, other: &Self) -> Ordering {
        if self.0 == other.0 {
            Ordering::Equal
        } else {
            self.0.total_cmp(&other.0)
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum GroupKey {
    Xic(Key),
    Srm(Key, Key),
}
struct Group {
    first: usize,
    points: Vec<ChromatogramPeak>,
}
struct Work {
    work: usize,
    bytes: usize,
}
impl Work {
    fn new(l: ChromatogramConversionLimits) -> Self {
        Self {
            work: l.max_work,
            bytes: l.max_bytes,
        }
    }
    fn meter(&mut self) -> Meter<'_> {
        Meter {
            work: &mut self.work,
            bytes: &mut self.bytes,
        }
    }
    fn spend(&mut self, n: usize) -> Result<()> {
        self.work = self.work.checked_sub(n).ok_or_else(limit)?;
        Ok(())
    }
    fn vector<T>(&mut self, n: usize) -> Result<Vec<T>> {
        self.bytes = self
            .bytes
            .checked_sub(mul(n, size_of::<T>())?)
            .ok_or_else(limit)?;
        let mut v = Vec::new();
        v.try_reserve_exact(n).map_err(|_| limit())?;
        Ok(v)
    }
    fn push<T>(&mut self, v: &mut Vec<T>, value: T) -> Result<()> {
        if v.len() == v.capacity() {
            let capacity = add(v.len(), 1)?.max(v.capacity().checked_mul(2).ok_or_else(limit)?);
            self.bytes = self
                .bytes
                .checked_sub(mul(capacity, size_of::<T>())?)
                .ok_or_else(limit)?;
            self.spend(v.len())?;
            v.try_reserve_exact(capacity - v.len())
                .map_err(|_| limit())?;
        }
        self.spend(1)?;
        v.push(value);
        Ok(())
    }
}
fn add(a: usize, b: usize) -> Result<usize> {
    a.checked_add(b).ok_or_else(limit)
}
fn mul(a: usize, b: usize) -> Result<usize> {
    a.checked_mul(b).ok_or_else(limit)
}
fn cap(n: usize, max: usize) -> Result<()> {
    if n > max { Err(limit()) } else { Ok(()) }
}
fn finite(v: f64) -> Result<()> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(Error::InvalidValue(
            "nonfinite chromatogram conversion coordinate/intensity".into(),
        ))
    }
}
fn limit() -> Error {
    Error::InvalidValue("chromatogram conversion resource limit exceeded".into())
}
