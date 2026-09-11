// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Format-independent streaming interfaces.

use crate::{MSChromatogram, MSSpectrum, Result, metadata::ExperimentalSettings};
use std::ops::ControlFlow;

/// Complete four-operation source IMSDataConsumer interface.
///
/// Records are borrowed mutably for one callback. A callback may change or move
/// their contents; no reference can outlive that call. Continue completes normal
/// consumption. Break requests a successful soft stop before the current record
/// is retained by its producer. An error stops production; earlier external
/// callback effects cannot be rolled back. There are no silent default methods.
pub trait MSDataConsumer {
    fn set_expected_size(&mut self, spectra: usize, chromatograms: usize) -> Result<()>;
    fn set_experimental_settings(&mut self, settings: &ExperimentalSettings) -> Result<()>;
    fn consume_spectrum(&mut self, spectrum: &mut MSSpectrum) -> Result<ControlFlow<()>>;
    fn consume_chromatogram(
        &mut self,
        chromatogram: &mut MSChromatogram,
    ) -> Result<ControlFlow<()>>;
}
