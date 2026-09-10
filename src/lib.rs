// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Native Rust mass-spectrometry data structures, chemistry, processing and I/O.
//!
//! This crate ports selected OpenMS core behavior. See `docs/PORTING_STATUS.md`
//! for the implemented API and explicit differences from the C++ library.
//!
//! ```
//! use openms::{MSSpectrum, Peak1D};
//! use openms::processing::{Normalizer, SpectrumFilter};
//! use openms::chemistry::AASequence;
//!
//! let mut spectrum = MSSpectrum::from_peaks(vec![
//!     Peak1D::new(100.0, 10.0), Peak1D::new(200.0, 40.0),
//! ]);
//! Normalizer::default().filter_spectrum(&mut spectrum)?;
//! assert_eq!(spectrum.peaks[0].intensity, 0.25);
//! assert_eq!(spectrum.find_nearest(199.0)?, Some(1));
//! let peptide = AASequence::parse("DFPIANGER")?;
//! assert!((peptide.mz(2)? - 509.751259).abs() < 1e-5);
//! # Ok::<(), openms::Error>(())
//! ```

pub mod analysis;
pub mod chemistry;
pub mod comparison;
pub mod error;
pub mod format;
pub mod identification;
pub mod kernel;
pub mod metadata;
pub mod processing;

pub use error::{Error, Result};
pub use kernel::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D, Precursor};

/// Version of the upstream scientific SDK targeted by this native Rust port.
/// This is distinct from the Rust crate version and does not imply full API parity.
pub const CORE_SDK_VERSION: &str = "4.0.0";
/// Exact upstream SDK source used for the current compatibility target.
/// Historical scientific fixtures retain their original, independently recorded pins.
pub const CORE_SDK_REVISION: &str = "6bfc0e4711105f4eda2fea86812a83af7c7e791f";
