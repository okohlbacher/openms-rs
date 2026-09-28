// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A consumer of spectra and chromatograms that performs no operation.
//!
//! Ports `FORMAT/DATAACCESS/NoopMSDataConsumer.h` (header-only apart from one
//! file-scope object in its `.cpp`). See `docs/MS_DATA_CONSUMERS_SUPPORT.md`
//! for the API mapping and the evidence.

use crate::interfaces::MSDataConsumer;
use crate::kernel::{MSChromatogram, MSSpectrum};
use crate::metadata::ExperimentalSettings;
use crate::Result;
use std::ops::ControlFlow;

/// Consumer that performs no operation.
///
/// This is sometimes necessary to fulfil the requirement of passing a valid
/// [`MSDataConsumer`] when no operation is required. Every call succeeds,
/// leaves its argument unchanged and asks the producer to continue.
///
/// The source's `NoopMSDataConsumer.cpp` also defines a file-scope object
/// `default_noop_consumer` that no header declares, so no caller can name it;
/// it has no counterpart here. A value of this zero-sized type costs nothing
/// to create where one is needed.
///
/// This type discards records without counting them. The streaming mzML
/// writer's own do-nothing variant, which also reports what arrived, is
/// [`NoopMSDataWritingConsumer`](crate::format::ms_data_writing_consumer::NoopMSDataWritingConsumer)
/// behind the `mzml` feature.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct NoopMSDataConsumer;

impl NoopMSDataConsumer {
    /// A consumer that does nothing; source `NoopMSDataConsumer()`.
    pub fn new() -> Self {
        Self
    }
}

impl MSDataConsumer for NoopMSDataConsumer {
    /// Ignores the expected sizes; source `setExpectedSize(Size, Size) {}`.
    fn set_expected_size(&mut self, _spectra: usize, _chromatograms: usize) -> Result<()> {
        Ok(())
    }

    /// Ignores the settings; source `setExperimentalSettings(const ExperimentalSettings&) {}`.
    fn set_experimental_settings(&mut self, _settings: &ExperimentalSettings) -> Result<()> {
        Ok(())
    }

    /// Leaves the spectrum unchanged and continues; source `consumeSpectrum(SpectrumType&) {}`.
    fn consume_spectrum(&mut self, _spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>> {
        Ok(ControlFlow::Continue(()))
    }

    /// Leaves the chromatogram unchanged and continues; source
    /// `consumeChromatogram(ChromatogramType&) {}`.
    fn consume_chromatogram(
        &mut self,
        _chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>> {
        Ok(ControlFlow::Continue(()))
    }
}
