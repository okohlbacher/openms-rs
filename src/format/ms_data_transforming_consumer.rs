// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that applies caller-supplied functions to MS data on the fly.
//!
//! Ports `FORMAT/DATAACCESS/MSDataTransformingConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.

use crate::Result;
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use std::fmt;
use std::ops::ControlFlow;

/// A function applied to every consumed spectrum.
pub type SpectrumFunc<'a> = Box<dyn FnMut(&mut MSSpectrum) -> Result<()> + 'a>;
/// A function applied to every consumed chromatogram.
pub type ChromatogramFunc<'a> = Box<dyn FnMut(&mut MSChromatogram) -> Result<()> + 'a>;
/// A function applied to the experimental settings when they are set.
pub type ExperimentalSettingsFunc<'a> = Box<dyn FnMut(&ExperimentalSettings) -> Result<()> + 'a>;

/// Transforming consumer of MS data.
///
/// Works on each spectrum or chromatogram on the fly, while it is being read,
/// using a function provided by the caller. The function can transform the
/// record in place, or just extract some information and store it in captured
/// variables.
///
/// The source requires the caller to "make sure the captured variables are
/// still in scope when this consumer is used"; here the lifetime `'a` of the
/// stored closures makes that a compile-time check, so a closure may borrow
/// local state mutably and the borrow ends when the consumer is dropped.
///
/// The source's functions return `void`; these return [`Result`], so a
/// function can fail the consumption instead of throwing. A function that
/// cannot fail returns `Ok(())`.
///
/// ```
/// use openms::format::ms_data_transforming_consumer::MSDataTransformingConsumer;
/// use openms::interfaces::MSDataConsumer;
/// use openms::kernel::{MSSpectrum, Peak1D};
///
/// let mut total_peaks = 0;
/// {
///     let mut consumer = MSDataTransformingConsumer::new();
///     consumer.set_spectra_processing_func(|s: &mut MSSpectrum| {
///         total_peaks += s.peaks.len();
///         Ok(())
///     });
///     let mut spectrum = MSSpectrum::default();
///     spectrum.peaks = vec![Peak1D::new(100.0, 1.0), Peak1D::new(200.0, 2.0)];
///     consumer.consume_spectrum(&mut spectrum)?;
/// }
/// assert_eq!(total_peaks, 2);
/// # Ok::<(), openms::Error>(())
/// ```
#[derive(Default)]
pub struct MSDataTransformingConsumer<'a> {
    /// Source `lambda_spec_`; `None` is the source's `nullptr`.
    spectrum_func: Option<SpectrumFunc<'a>>,
    /// Source `lambda_chrom_`.
    chromatogram_func: Option<ChromatogramFunc<'a>>,
    /// Source `lambda_exp_settings_`.
    settings_func: Option<ExperimentalSettingsFunc<'a>>,
}

impl fmt::Debug for MSDataTransformingConsumer<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("MSDataTransformingConsumer")
            .field("spectrum_func", &self.spectrum_func.is_some())
            .field("chromatogram_func", &self.chromatogram_func.is_some())
            .field("settings_func", &self.settings_func.is_some())
            .finish()
    }
}

impl<'a> MSDataTransformingConsumer<'a> {
    /// A consumer with no function set, which leaves every record unchanged;
    /// source `MSDataTransformingConsumer()`, whose three members start as
    /// `nullptr`.
    pub fn new() -> Self {
        Self::default()
    }

    /// Sets the function called for every spectrum passed to this consumer;
    /// source `setSpectraProcessingFunc`.
    ///
    /// The function may capture local state, which lets it change that state;
    /// the borrow lasts as long as the consumer. A later call replaces the
    /// function.
    pub fn set_spectra_processing_func(
        &mut self,
        func: impl FnMut(&mut MSSpectrum) -> Result<()> + 'a,
    ) {
        self.spectrum_func = Some(Box::new(func));
    }

    /// Removes the spectrum function, so spectra are left unchanged; source
    /// `setSpectraProcessingFunc(nullptr)`.
    pub fn clear_spectra_processing_func(&mut self) {
        self.spectrum_func = None;
    }

    /// Sets the function called for every chromatogram passed to this
    /// consumer; source `setChromatogramProcessingFunc`.
    ///
    /// The same capture rules apply as for
    /// [`set_spectra_processing_func`](Self::set_spectra_processing_func).
    pub fn set_chromatogram_processing_func(
        &mut self,
        func: impl FnMut(&mut MSChromatogram) -> Result<()> + 'a,
    ) {
        self.chromatogram_func = Some(Box::new(func));
    }

    /// Removes the chromatogram function, so chromatograms are left
    /// unchanged; source `setChromatogramProcessingFunc(nullptr)`.
    pub fn clear_chromatogram_processing_func(&mut self) {
        self.chromatogram_func = None;
    }

    /// Sets the function called when the experimental settings are set through
    /// this consumer; source `setExperimentalSettingsFunc`.
    ///
    /// Without one (the default) setting the settings does nothing.
    pub fn set_experimental_settings_func(
        &mut self,
        func: impl FnMut(&ExperimentalSettings) -> Result<()> + 'a,
    ) {
        self.settings_func = Some(Box::new(func));
    }

    /// Removes the settings function; source
    /// `setExperimentalSettingsFunc(nullptr)`.
    pub fn clear_experimental_settings_func(&mut self) {
        self.settings_func = None;
    }

    /// Whether a spectrum function is set.
    ///
    /// No source counterpart; the source's `lambda_spec_` is protected.
    pub fn has_spectra_processing_func(&self) -> bool {
        self.spectrum_func.is_some()
    }

    /// Whether a chromatogram function is set; protected `lambda_chrom_` upstream.
    pub fn has_chromatogram_processing_func(&self) -> bool {
        self.chromatogram_func.is_some()
    }

    /// Whether a settings function is set; protected `lambda_exp_settings_` upstream.
    pub fn has_experimental_settings_func(&self) -> bool {
        self.settings_func.is_some()
    }
}

impl MSDataConsumer for MSDataTransformingConsumer<'_> {
    /// Ignored, as in the source.
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }

    /// Passes the settings to the settings function, if one is set; source
    /// `setExperimentalSettings`.
    ///
    /// The header documents this override as "ignored", but the implementation
    /// calls `lambda_exp_settings_` when it is set; this follows the
    /// implementation.
    ///
    /// # Errors
    ///
    /// The function's error.
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        match self.settings_func.as_mut() {
            Some(func) => func(settings),
            None => Ok(()),
        }
    }

    /// Applies the spectrum function to the spectrum in place, if one is set,
    /// and continues; source `consumeSpectrum`.
    ///
    /// # Errors
    ///
    /// The function's error.
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        if let Some(func) = self.spectrum_func.as_mut() {
            func(spectrum)?;
        }
        Ok(ControlFlow::Continue(()))
    }

    /// Applies the chromatogram function to the chromatogram in place, if one
    /// is set, and continues; source `consumeChromatogram`.
    ///
    /// # Errors
    ///
    /// The function's error.
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        if let Some(func) = self.chromatogram_func.as_mut() {
            func(chromatogram)?;
        }
        Ok(ControlFlow::Continue(()))
    }
}
