// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that merges spectra sharing a retention time before passing them on.
//!
//! Ports `FORMAT/DATAACCESS/MSDataAggregatingConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.
//!
//! The source merges each group with `SpectrumAddition::addUpSpectra(list, -1,
//! true)` from `ANALYSIS/OPENSWATH/SpectrumAddition.h`. That header is not
//! ported yet, and `format` may not depend on `analysis` in any case, so the
//! merge is a function the caller supplies; see
//! [`MSDataAggregatingConsumer::new`](crate::format::ms_data_aggregating_consumer::MSDataAggregatingConsumer::new).

use crate::Result;
use crate::interfaces::MSDataConsumer;
use crate::kernel::spectrum_helper::copy_spectrum_meta;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use std::fmt;
use std::ops::ControlFlow;

/// Two consecutive spectra whose retention times differ by less than this
/// (seconds, strictly) belong to the same group.
///
/// Source `std::fabs(RT - previous_rt_) < 1e-5` in `consumeSpectrum`. A NaN
/// retention time compares false, so it always starts a new group, as it does
/// in the source.
pub const RT_TOLERANCE: f64 = 1e-5;

/// Aggregates spectra by retention time.
///
/// Consecutive spectra passed to this consumer that have the same retention
/// time (within [`RT_TOLERANCE`] of the previous spectrum's) are merged, and
/// the merged spectrum is passed to the next consumer. Chromatograms are passed
/// on unchanged, immediately.
///
/// A group is passed on when the first spectrum of the next group arrives, or
/// by [`finish`](Self::finish) for the last group. Each merged spectrum
/// carries the metadata of the group's *first* spectrum (source
/// `copySpectrumMeta(s_list[0], tmps, false)`), with the peaks and data arrays
/// the merge function produced. A group of one spectrum is merged too, exactly
/// as in the source, so the merge function decides whether a single spectrum
/// passes through unchanged.
///
/// The next consumer is borrowed, not owned, as the source's constructor
/// "does not transfer ownership of the consumer"; the source's warning that
/// deleting the next consumer before this one risks a memory error is what the
/// lifetime `'a` rules out.
///
/// The source does **not** forward `setExpectedSize` or
/// `setExperimentalSettings` to the next consumer: both are empty overrides.
/// This port reproduces that, so a caller that needs the next consumer
/// configured configures it directly before chaining it here.
pub struct MSDataAggregatingConsumer<'a, F>
where
    F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>,
{
    /// Source `next_consumer_`.
    next_consumer: &'a mut dyn MSDataConsumer,
    /// The merge the source hard-codes as `SpectrumAddition::addUpSpectra(list, -1, true)`.
    add_up: F,
    /// Source `previous_rt_`.
    previous_rt: f64,
    /// Source `rt_initialized_`.
    rt_initialized: bool,
    /// The spectra of the current group; source `s_list`.
    pending: Vec<MSSpectrum>,
    /// Set when the next consumer asked to stop; nothing more is passed on.
    stopped: bool,
}

impl<F> fmt::Debug for MSDataAggregatingConsumer<'_, F>
where
    F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MSDataAggregatingConsumer")
            .field("previous_rt", &self.previous_rt)
            .field("rt_initialized", &self.rt_initialized)
            .field("pending", &self.pending.len())
            .field("stopped", &self.stopped)
            .finish()
    }
}

impl<'a, F> MSDataAggregatingConsumer<'a, F>
where
    F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>,
{
    /// A consumer passing merged spectra to `next_consumer`; source
    /// `MSDataAggregatingConsumer(IMSDataConsumer* next_consumer)`.
    ///
    /// `add_up` merges one group of spectra (at least one, all with the same
    /// retention time) into one spectrum; only its peaks and data arrays are
    /// used, the metadata is copied from the group's first spectrum afterwards.
    /// The source's merge is `SpectrumAddition::addUpSpectra(spectra, -1,
    /// true)` — resampling to the auto-detected sampling rate and dropping
    /// zero-intensity points — which is not ported (the header is
    /// `ANALYSIS/OPENSWATH/SpectrumAddition.h`); until it is, the caller
    /// supplies it or another merge.
    ///
    /// The source's class test constructs one with a null next consumer
    /// ("don't do that ..."), which crashes on the first flush; a reference
    /// cannot be null.
    pub fn new(next_consumer: &'a mut dyn MSDataConsumer, add_up: F) -> Self {
        Self {
            next_consumer,
            add_up,
            previous_rt: 0.0,
            rt_initialized: false,
            pending: Vec::new(),
            stopped: false,
        }
    }

    /// Number of spectra waiting in the current, not yet passed-on group.
    ///
    /// No source counterpart; the source's `s_list` is private.
    pub fn pending(&self) -> usize {
        self.pending.len()
    }

    /// Passes the last group on and ends the consumption; the flush of the
    /// source's destructor, which "flushes data to next consumer".
    ///
    /// Returns the next consumer's answer for the last merged spectrum, or
    /// `Continue` when nothing was pending.
    ///
    /// Dropping the consumer without calling this performs the same flush,
    /// but a drop cannot report an error, so it logs a warning instead; call
    /// this to see the outcome.
    ///
    /// # Errors
    ///
    /// The merge function's error, or the next consumer's. The group is then
    /// discarded, and the drop that follows does not try again.
    pub fn finish(mut self) -> Result<ControlFlow<()>> {
        let result = self.flush();
        self.pending.clear();
        result
    }

    /// Merges and passes on the pending group, leaving it empty.
    fn flush(&mut self) -> Result<ControlFlow<()>> {
        if self.pending.is_empty() || self.stopped {
            self.pending.clear();
            return Ok(ControlFlow::Continue(()));
        }
        let mut merged = (self.add_up)(&self.pending)?;
        copy_spectrum_meta(&self.pending[0], &mut merged, false);
        let flow = self.next_consumer.consume_spectrum(&mut merged)?;
        self.pending.clear();
        if flow.is_break() {
            self.stopped = true;
        }
        Ok(flow)
    }
}

impl<F> Drop for MSDataAggregatingConsumer<'_, F>
where
    F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>,
{
    fn drop(&mut self) {
        if let Err(error) = self.flush() {
            // The source's destructor flushes too, and has nowhere to report
            // a failure either; a drop must not panic.
            let _ = crate::openms_log_warn!(
                "MSDataAggregatingConsumer: the last spectrum group was not passed on: {error}"
            );
        }
    }
}

impl<F> MSDataConsumer for MSDataAggregatingConsumer<'_, F>
where
    F: FnMut(&[MSSpectrum]) -> Result<MSSpectrum>,
{
    /// Ignored, and not forwarded to the next consumer, as in the source.
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }

    /// Ignored, and not forwarded to the next consumer, as in the source.
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }

    /// Adds the spectrum to the current group, or passes the current group on
    /// and starts a new one; source `consumeSpectrum`.
    ///
    /// The spectrum joins the current group when a spectrum was seen before
    /// and its retention time is within [`RT_TOLERANCE`] of the *previous*
    /// spectrum's (not the group's first), so a slow drift can chain a long
    /// group, as in the source. A copy is kept; the caller's spectrum is not
    /// changed.
    ///
    /// Returns the next consumer's answer when a group was passed on, and
    /// `Continue` otherwise. After the next consumer answered
    /// [`ControlFlow::Break`], this spectrum is not kept and every later one
    /// is refused the same way, since the producer was asked to stop.
    ///
    /// # Errors
    ///
    /// The merge function's error or the next consumer's, when a group is
    /// passed on. The group, this spectrum and the retention-time state are
    /// then as they were before the call, as when the source's call throws.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        if self.stopped {
            return Ok(ControlFlow::Break(()));
        }
        let rt = spectrum.rt;
        if self.rt_initialized && (rt - self.previous_rt).abs() < RT_TOLERANCE {
            self.pending.push(spectrum.clone());
        } else {
            if self.rt_initialized && !self.pending.is_empty() {
                let mut merged = (self.add_up)(&self.pending)?;
                copy_spectrum_meta(&self.pending[0], &mut merged, false);
                let flow = self.next_consumer.consume_spectrum(&mut merged)?;
                self.pending.clear();
                if flow.is_break() {
                    self.stopped = true;
                    return Ok(ControlFlow::Break(()));
                }
            }
            self.pending.clear();
            self.pending.push(spectrum.clone());
        }
        self.previous_rt = rt;
        self.rt_initialized = true;
        Ok(ControlFlow::Continue(()))
    }

    /// Passes the chromatogram to the next consumer unchanged; source
    /// `consumeChromatogram`.
    ///
    /// # Errors
    ///
    /// The next consumer's error.
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        if self.stopped {
            return Ok(ControlFlow::Break(()));
        }
        let flow = self.next_consumer.consume_chromatogram(chromatogram)?;
        if flow.is_break() {
            self.stopped = true;
        }
        Ok(flow)
    }
}
