// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer that keeps every spectrum and chromatogram passed to it in memory.
//!
//! Ports `FORMAT/DATAACCESS/MSDataStoringConsumer.h` and its `.cpp`. See
//! `docs/MS_DATA_CONSUMERS_SUPPORT.md` for the API mapping and the evidence.

use crate::Result;
use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use std::ops::ControlFlow;

/// The most records one [`MSDataConsumer::set_expected_size`] call reserves
/// room for, per list.
///
/// The expected size is a capacity hint: the source forwards it to
/// `reserveSpaceSpectra`/`reserveSpaceChromatograms`, which only reserve. A
/// producer passes whatever count its input declares, so an unchecked
/// reservation would let one attribute of an untrusted file request an
/// arbitrary allocation before a single record exists. Reserving at most this
/// many slots keeps the hint's benefit for ordinary runs; a larger run still
/// stores every record, the vector simply grows as it is filled.
pub const MAX_RESERVED_RECORDS: usize = 1 << 16;

/// Consumer that simply stores the data.
///
/// Keeps the spectra and chromatograms passed to it in memory, in arrival
/// order; the data is accessed through [`MSDataStoringConsumer::data`].
///
/// Each record is stored as a copy of its state when it arrived, as the source
/// copies it into its `PeakMap` (`exp_.addSpectrum(s)`); the caller's record is
/// left as it was.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MSDataStoringConsumer {
    /// The stored experiment; source `exp_`.
    experiment: MSExperiment,
}

impl MSDataStoringConsumer {
    /// An empty consumer; source `MSDataStoringConsumer()`.
    pub fn new() -> Self {
        Self::default()
    }

    /// The data stored so far; source `getData()`.
    pub fn data(&self) -> &MSExperiment {
        &self.experiment
    }

    /// Consumes the consumer and returns the stored data without copying it.
    ///
    /// No source counterpart: the source keeps its `PeakMap` inside the
    /// consumer and hands out a `const&`.
    pub fn into_data(self) -> MSExperiment {
        self.experiment
    }
}

impl MSDataConsumer for MSDataStoringConsumer {
    /// Reserves room for the expected records; source `setExpectedSize`.
    ///
    /// Each list reserves at most [`MAX_RESERVED_RECORDS`] slots, and a failed
    /// reservation is ignored, because the source's call only reserves and a
    /// reservation changes no stored value. Always succeeds.
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()> {
        let spectra = spectra.min(MAX_RESERVED_RECORDS);
        let chromatograms = chromatograms.min(MAX_RESERVED_RECORDS);
        // A reservation failure changes nothing observable; the vectors grow on
        // demand instead, as they do past the cap.
        let _ = self.experiment.spectra.try_reserve(spectra);
        let _ = self.experiment.chromatograms.try_reserve(chromatograms);
        Ok(())
    }

    /// Replaces the stored experiment's settings and keeps its data.
    ///
    /// Source `exp_ = settings`, which assigns only the `ExperimentalSettings`
    /// base of the `PeakMap` ("only override the settings, keep the data"), so
    /// records consumed before or after the call are all kept. A later call
    /// replaces the settings again.
    ///
    /// The source's `sqMassRunID` travels as a meta value of the settings; here
    /// it is [`MSExperiment::sql_run_id`], outside
    /// [`ExperimentalSettings`], so this call does not touch it and it stays
    /// zero, the value the source reports when the settings carry none.
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()> {
        self.experiment.settings = settings.clone();
        Ok(())
    }

    /// Appends a copy of the spectrum; source `consumeSpectrum` (`exp_.addSpectrum(s)`).
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        self.experiment.spectra.push(spectrum.clone());
        Ok(ControlFlow::Continue(()))
    }

    /// Appends a copy of the chromatogram; source `consumeChromatogram`
    /// (`exp_.addChromatogram(c)`).
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        self.experiment.chromatograms.push(chromatogram.clone());
        Ok(ControlFlow::Continue(()))
    }
}
