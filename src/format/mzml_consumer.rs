// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Two-pass setup and bounded, independently flushed record pools.

use super::*;
use crate::metadata::{ExperimentalSettings, MetaInfo};
use std::{mem::size_of, ops::ControlFlow, path::Path};

/// Scientific choices plus administrative allocation/work limits.
/// Parser, binary and selection limits apply to each pass as documented by
/// ReadOptions/LoadOptions. They never reset at a pool boundary. Callback work
/// and callback-created payload are caller-owned, not charged by this interface.
#[derive(Clone, Debug)]
pub struct TransformOptions {
    pub skip_full_count: bool,
    pub skip_first_pass: bool,
    pub load: LoadOptions,
    pub read: ReadOptions,
    /// Cumulative own settings-copy, pool-vector and retaining-vector storage.
    pub max_bytes: usize,
    /// Cumulative own settings-copy and descriptor-move work.
    /// Source-file deduplication uses the separate per-pass header allowance.
    pub max_work: usize,
}
impl Default for TransformOptions {
    fn default() -> Self {
        Self {
            skip_full_count: false,
            skip_first_pass: false,
            load: LoadOptions::default(),
            read: ReadOptions::default(),
            max_bytes: 256 * 1024 * 1024,
            max_work: 50_000_000,
        }
    }
}
/// Report for successful completion or a callback-requested soft stop.
/// An error may follow earlier callbacks; their effects cannot be rolled back.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct TransformReport {
    /// None when the setup pass was suppressed; raw declarations otherwise.
    pub expected: Option<MzMLCounts>,
    /// Callbacks completed, including the callback requesting a soft stop.
    pub delivered: MzMLCounts,
    pub stopped: bool,
}

pub fn transform(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
) -> Result<TransformReport> {
    transform_with_options(path, consumer, &TransformOptions::default())
}
pub fn transform_with_options(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
    options: &TransformOptions,
) -> Result<TransformReport> {
    let path = path.as_ref();
    transform_from(|| crate::format::path_io::open(path), consumer, options)
}
pub fn transform_into(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
    destination: &mut MSExperiment,
) -> Result<TransformReport> {
    transform_into_with_options(path, consumer, destination, &TransformOptions::default())
}
/// Append callback-modified records only after successful parsing or soft stop.
/// Prior record payloads are moved only by vector growth, never deeply cloned.
/// Error leaves logical destination contents intact, but callbacks and a final
/// reservation's capacity changes are external observable effects.
pub fn transform_into_with_options(
    path: impl AsRef<Path>,
    consumer: &mut dyn MSDataConsumer,
    destination: &mut MSExperiment,
    options: &TransformOptions,
) -> Result<TransformReport> {
    let path = path.as_ref();
    transform_from_into(
        || crate::format::path_io::open(path),
        consumer,
        destination,
        options,
    )
}
/// The factory is called once per parser pass, twice normally or once when
/// skip_first_pass is true. It permits compressed/nonseekable inputs without a
/// hidden whole-document buffer. Distinct results between opens remain the
/// caller's responsibility, as with files changed between source passes.
pub fn transform_from<R: BufRead>(
    open: impl FnMut() -> Result<R>,
    consumer: &mut dyn MSDataConsumer,
    options: &TransformOptions,
) -> Result<TransformReport> {
    run(open, consumer, None, options)
}
pub fn transform_from_into<R: BufRead>(
    open: impl FnMut() -> Result<R>,
    consumer: &mut dyn MSDataConsumer,
    destination: &mut MSExperiment,
    options: &TransformOptions,
) -> Result<TransformReport> {
    run(open, consumer, Some(destination), options)
}
fn resource() -> Error {
    Error::InvalidValue("mzML consumer resource limit exceeded".into())
}
pub(super) struct Admin {
    work: usize,
    bytes: usize,
}
impl Admin {
    fn spend(&mut self, work: usize, bytes: usize) -> Result<()> {
        self.work = self.work.checked_sub(work).ok_or_else(resource)?;
        self.bytes = self.bytes.checked_sub(bytes).ok_or_else(resource)?;
        Ok(())
    }
    fn slots<T>(&mut self, n: usize) -> Result<()> {
        self.spend(n, n.checked_mul(size_of::<T>()).ok_or_else(resource)?)
    }
    fn push<T>(&mut self, values: &mut Vec<T>, value: T, maximum: usize) -> Result<()> {
        self.spend(1, 0)?;
        if values.len() == values.capacity() {
            let capacity = values
                .capacity()
                .checked_mul(2)
                .ok_or_else(resource)?
                .max(4)
                .min(maximum.max(1));
            if capacity <= values.len() {
                return Err(resource());
            }
            self.slots::<T>(capacity)?;
            self.spend(values.len(), 0)?;
            values
                .try_reserve_exact(capacity - values.len())
                .map_err(|_| resource())?;
        }
        values.push(value);
        Ok(())
    }
    fn reserve_append<T>(&mut self, old: &mut Vec<T>, added: usize) -> Result<()> {
        let total = old.len().checked_add(added).ok_or_else(resource)?;
        self.spend(added, 0)?;
        if total > old.capacity() {
            self.slots::<T>(total)?;
            self.spend(old.len(), 0)?;
            old.try_reserve_exact(added).map_err(|_| resource())?;
        }
        Ok(())
    }
}
fn run<R: BufRead>(
    mut open: impl FnMut() -> Result<R>,
    consumer: &mut dyn MSDataConsumer,
    destination: Option<&mut MSExperiment>,
    options: &TransformOptions,
) -> Result<TransformReport> {
    options.load.validate()?;
    let mut admin = Admin {
        work: options.max_work,
        bytes: options.max_bytes,
    };
    admin.spend(1, 0)?;
    // Validate/copy only the settings that will participate in source merging.
    // Old spectra/chromatograms and their arbitrary graphs are not traversed.
    let mut settings = if let Some(old) = &destination {
        old.settings
            .with_budget(&mut admin.work, &mut admin.bytes)?;
        old.settings.clone()
    } else {
        ExperimentalSettings::default()
    };
    let mut previous_metadata = std::mem::take(&mut settings.metadata);
    let mut report = TransformReport::default();
    if !options.skip_first_pass {
        let (settings, counts) = counts::setup(
            open()?,
            &options.load.scientific,
            &options.read,
            options.skip_full_count,
        )?;
        consumer.set_expected_size(counts.spectra, counts.chromatograms)?;
        consumer.set_experimental_settings(&settings)?;
        report.expected = Some(counts);
    }
    let retaining = destination.is_some();
    let mut sink = Sink {
        consumer,
        admin: &mut admin,
        report: &mut report,
        maximum: options.load.scientific.max_data_pool_size.max(1),
        record_limit: options.read.max_records,
        spectra: Vec::new(),
        chromatograms: Vec::new(),
        retained: retaining.then(|| (Vec::new(), Vec::new())),
    };
    let parsed = super::read_engine(
        open()?,
        &options.read,
        Some(&options.load),
        options.load.scientific.metadata_only,
        MSExperiment {
            settings,
            ..Default::default()
        },
        Some(&mut sink),
    )?;
    let retained = sink.retained.take();
    drop(sink);
    if let Some(destination) = destination {
        let mut settings = parsed.settings;
        // Old keys are overwritten by present input keys; duplicate input keys
        // still pass through the reader's ordinary duplicate checks.
        merge_metadata(&mut previous_metadata, &mut settings.metadata, &mut admin)?;
        settings.metadata = previous_metadata;
        let (mut spectra, mut chromatograms) = retained.unwrap();
        admin.reserve_append(&mut destination.spectra, spectra.len())?;
        admin.reserve_append(&mut destination.chromatograms, chromatograms.len())?;
        destination.spectra.append(&mut spectra);
        destination.chromatograms.append(&mut chromatograms);
        destination.settings = settings;
    }
    Ok(report)
}
fn merge_metadata(old: &mut MetaInfo, new: &mut MetaInfo, admin: &mut Admin) -> Result<()> {
    let n = old.len().checked_add(new.len()).ok_or_else(resource)?;
    if n == 0 {
        return Ok(());
    }
    admin.spend(n, 0)?;
    let width = old
        .keys()
        .chain(new.keys())
        .map(String::len)
        .max()
        .unwrap_or(0);
    let comparisons = (usize::BITS - n.leading_zeros()) as usize;
    admin.spend(
        n.saturating_mul(comparisons)
            .saturating_mul(12)
            .saturating_mul(width.saturating_add(1)),
        0,
    )?;
    crate::kernel::data_array::Meter {
        work: &mut admin.work,
        bytes: &mut admin.bytes,
    }
    .tree::<(String, crate::metadata::MetaValue)>(n)?;
    old.append(new);
    Ok(())
}

// This is one transient stack value, immediately moved into its record Vec.
// Boxing would add an otherwise unnecessary allocation for every chromatogram.
#[allow(clippy::large_enum_variant)]
pub(super) enum Completed {
    Spectrum(MSSpectrum),
    Chromatogram(MSChromatogram),
}
pub(super) struct Sink<'a> {
    consumer: &'a mut dyn MSDataConsumer,
    admin: &'a mut Admin,
    report: &'a mut TransformReport,
    maximum: usize,
    record_limit: usize,
    spectra: Vec<MSSpectrum>,
    chromatograms: Vec<MSChromatogram>,
    retained: Option<(Vec<MSSpectrum>, Vec<MSChromatogram>)>,
}
impl Sink<'_> {
    pub(super) fn accept(&mut self, record: Completed) -> Result<bool> {
        match record {
            Completed::Spectrum(s) => {
                self.admin
                    .push(&mut self.spectra, s, self.maximum.min(self.record_limit))?;
                if self.spectra.len() >= self.maximum {
                    self.flush_spectra()
                } else {
                    Ok(true)
                }
            }
            Completed::Chromatogram(c) => {
                self.admin.push(
                    &mut self.chromatograms,
                    c,
                    self.maximum.min(self.record_limit),
                )?;
                if self.chromatograms.len() >= self.maximum {
                    self.flush_chromatograms()
                } else {
                    Ok(true)
                }
            }
        }
    }
    fn flush_spectra(&mut self) -> Result<bool> {
        for mut spectrum in self.spectra.drain(..) {
            self.admin.spend(1, 0)?;
            let control = self.consumer.consume_spectrum(&mut spectrum)?;
            self.report.delivered.spectra += 1; // bounded by max_records
            if matches!(control, ControlFlow::Break(())) {
                self.report.stopped = true;
                return Ok(false);
            }
            if let Some((spectra, _)) = &mut self.retained {
                self.admin.push(spectra, spectrum, self.record_limit)?;
            }
        }
        Ok(true)
    }
    fn flush_chromatograms(&mut self) -> Result<bool> {
        for mut chromatogram in self.chromatograms.drain(..) {
            self.admin.spend(1, 0)?;
            let control = self.consumer.consume_chromatogram(&mut chromatogram)?;
            self.report.delivered.chromatograms += 1;
            if matches!(control, ControlFlow::Break(())) {
                self.report.stopped = true;
                return Ok(false);
            }
            if let Some((_, chromatograms)) = &mut self.retained {
                self.admin
                    .push(chromatograms, chromatogram, self.record_limit)?;
            }
        }
        Ok(true)
    }
    pub(super) fn finish(&mut self) -> Result<bool> {
        if !self.flush_spectra()? {
            return Ok(false);
        }
        self.flush_chromatograms()
    }
}
