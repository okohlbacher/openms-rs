// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that passes all consumed data through a sequence of consumers.
//!
//! Ports `FORMAT/DATAACCESS/MSDataChainingConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.

use crate::Result;
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use std::fmt;
use std::ops::ControlFlow;

/// Consumer that passes all consumed data through a set of operations.
///
/// Chains several consumers and applies them in a pre-specified order. This
/// is useful when an operation on a data set needs some pre-processing (data
/// reduction and the like) or post-processing (writing to disk, caching on
/// disk): the steps are added to the chain in the correct order, and the chain
/// is passed wherever an [`MSDataConsumer`] is expected, without that code
/// knowing about the individual steps.
///
/// The chain borrows its consumers and does not own them, as the source's
/// holds raw pointers it never deletes. The borrow ends when the chain is
/// dropped, after which each consumer is usable on its own again; the source's
/// "the caller must keep every consumer alive and delete it afterwards" is
/// what the lifetime `'a` checks at compile time.
///
/// ```
/// use openms::format::ms_data_chaining_consumer::MSDataChainingConsumer;
/// use openms::format::ms_data_storing_consumer::MSDataStoringConsumer;
/// use openms::format::ms_data_transforming_consumer::MSDataTransformingConsumer;
/// use openms::interfaces::MSDataConsumer;
/// use openms::kernel::{MSSpectrum, Peak1D};
///
/// let mut sort = MSDataTransformingConsumer::new();
/// sort.set_spectra_processing_func(|s: &mut MSSpectrum| s.sort_by_intensity(false));
/// let mut storage = MSDataStoringConsumer::new();
/// {
///     let steps: Vec<&mut dyn MSDataConsumer> = vec![&mut sort, &mut storage];
///     let mut chain = MSDataChainingConsumer::from_consumers(steps);
///     let mut spectrum = MSSpectrum::default();
///     spectrum.peaks = vec![Peak1D::new(1.0, 5.0), Peak1D::new(2.0, 1.0)];
///     chain.set_expected_size(1, 0)?;
///     chain.consume_spectrum(&mut spectrum)?;
/// }
/// // The storing consumer saw the spectrum after the sort.
/// assert_eq!(storage.data().spectra[0].peaks[0].mz, 2.0);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Default)]
pub struct MSDataChainingConsumer<'a> {
    /// The chained consumers, in call order; source `consumers_`.
    consumers: Vec<&'a mut dyn MSDataConsumer>,
}

impl fmt::Debug for MSDataChainingConsumer<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MSDataChainingConsumer")
            .field("consumers", &self.consumers.len())
            .finish()
    }
}

impl<'a> MSDataChainingConsumer<'a> {
    /// An empty chain; source `MSDataChainingConsumer()`.
    ///
    /// An empty chain accepts every call and does nothing.
    pub fn new() -> Self {
        Self::default()
    }

    /// A chain calling `consumers` sequentially, first to last; source
    /// `MSDataChainingConsumer(std::vector<IMSDataConsumer*>)`.
    ///
    /// Ownership is not transferred: the chain borrows each consumer.
    pub fn from_consumers(consumers: Vec<&'a mut dyn MSDataConsumer>) -> Self {
        Self { consumers }
    }

    /// Appends a consumer to the end of the chain; source `appendConsumer`.
    ///
    /// Ownership is not transferred: the chain borrows the consumer.
    pub fn append_consumer(&mut self, consumer: &'a mut dyn MSDataConsumer) {
        self.consumers.push(consumer);
    }

    /// Number of chained consumers.
    ///
    /// No source counterpart; the source's `consumers_` is private.
    pub fn len(&self) -> usize {
        self.consumers.len()
    }

    /// Whether the chain holds no consumer.
    pub fn is_empty(&self) -> bool {
        self.consumers.is_empty()
    }
}

impl MSDataConsumer for MSDataChainingConsumer<'_> {
    /// Sets the expected size for all chained consumers, in chain order;
    /// source `setExpectedSize`.
    ///
    /// # Errors
    ///
    /// The first consumer's error; the consumers after it are not called, as
    /// an exception leaves the source's loop.
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        for consumer in &mut self.consumers {
            consumer.set_expected_size(spectra, chromatograms)?;
        }
        Ok(())
    }

    /// Sets the experimental settings for all chained consumers, in chain
    /// order; source `setExperimentalSettings`.
    ///
    /// # Errors
    ///
    /// The first consumer's error; the consumers after it are not called.
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        for consumer in &mut self.consumers {
            consumer.set_experimental_settings(settings)?;
        }
        Ok(())
    }

    /// Calls all consumers in the specified order for the given spectrum;
    /// source `consumeSpectrum`.
    ///
    /// Every consumer receives the same spectrum, so a change one consumer
    /// makes is what the next one sees, and the caller sees the result of the
    /// whole chain.
    ///
    /// # Errors
    ///
    /// The first consumer's error; the consumers after it are not called.
    ///
    /// A consumer answering [`ControlFlow::Break`] asks its producer to stop
    /// before the record is retained. The source's interface has no soft stop,
    /// so this is native: the chain passes the record no further and answers
    /// `Break` itself, so no later consumer — a writer at the end of the chain,
    /// typically — retains a record an earlier step refused.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        for consumer in &mut self.consumers {
            if consumer.consume_spectrum(spectrum)?.is_break() {
                return Ok(ControlFlow::Break(()));
            }
        }
        Ok(ControlFlow::Continue(()))
    }

    /// Calls all consumers in the specified order for the given chromatogram;
    /// source `consumeChromatogram`.
    ///
    /// # Errors
    ///
    /// As [`consume_spectrum`](Self::consume_spectrum), including the native
    /// handling of [`ControlFlow::Break`].
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        for consumer in &mut self.consumers {
            if consumer.consume_chromatogram(chromatogram)?.is_break() {
                return Ok(ControlFlow::Break(()));
            }
        }
        Ok(ControlFlow::Continue(()))
    }
}
