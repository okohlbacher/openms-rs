// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Class-test port for `KERNEL/DPeak.h` and `KERNEL/StandardTypes.h`
//! (source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4, tier 3).
//!
//! Both headers declare names, not behaviour. `DPeak<1>::Type` / `DPeak<2>::Type`
//! are the concrete `Peak1D` / `Peak2D` in Rust, and the three `StandardTypes.h`
//! typedefs are `MSSpectrum`, `MSExperiment` and `MSChromatogram` themselves.
//! The source class tests only construct and destroy each name; every section
//! is reproduced here so the mapping is exercised rather than asserted.

use openms::comparison::{BinConfig, BinnedSpectrum};
use openms::kernel::Peak2D;
use openms::{ChromatogramPeak, MSChromatogram, MSExperiment, MSSpectrum, Peak1D};

/// `DPeak_test.cpp` sections `DPeak()` and `~DPeak()`: `DPeak<1>::Type` is
/// `Peak1D`; default construction succeeds and the value drops.
#[test]
fn dpeak_1_type_is_peak1d() {
    let peak: Peak1D = Peak1D::default();
    assert_eq!(peak, Peak1D::new(0.0, 0.0));
    assert_eq!(peak.mz, 0.0);
    assert_eq!(peak.intensity, 0.0);
    // Source `Peak1D::PositionType` is `DPosition<1>`; Rust stores the scalar.
    let boxed: Box<Peak1D> = Box::default();
    assert_eq!(*boxed, peak);
    drop(boxed);
}

/// `DPeak_test.cpp` sections `[EXTRA]DPeak()` and `[EXTRA]~DPeak()`:
/// `DPeak<2>::Type` is `Peak2D`.
#[test]
fn dpeak_2_type_is_peak2d() {
    let peak: Peak2D = Peak2D::default();
    assert_eq!(peak, Peak2D::new(0.0, 0.0, 0.0));
    assert_eq!(peak.position, [0.0, 0.0]);
    assert_eq!(Peak2D::DIMENSION, 2);
    let boxed: Box<Peak2D> = Box::default();
    assert_eq!(*boxed, peak);
    drop(boxed);
}

/// `DPeak.cpp` instantiates one global default `DPeak<1>`, `DPeak<2>` and one
/// of each `::Type`. Rust has no metafunction object to instantiate; the two
/// concrete defaults are the only observable state, and they are zero.
#[test]
fn dpeak_cpp_globals_have_only_the_two_concrete_defaults() {
    assert_eq!(
        Peak1D::default(),
        Peak1D {
            mz: 0.0,
            intensity: 0.0
        }
    );
    assert_eq!(
        Peak2D::default(),
        Peak2D {
            position: [0.0, 0.0],
            intensity: 0.0
        }
    );
}

/// `StandardTypes_test.cpp` `GOOD_TYPEDEF(PeakSpectrum)` (declared twice in the
/// source): `PeakSpectrum` is `MSSpectrum`. The source's `StandardTypes`
/// section is `NOT_TESTABLE`; the typedef check is construct-and-drop.
#[test]
fn peak_spectrum_typedef_is_ms_spectrum() {
    let spectrum: MSSpectrum = MSSpectrum::new();
    assert_eq!(spectrum, MSSpectrum::default());
    assert!(spectrum.peaks.is_empty());
    assert_eq!(spectrum.rt, -1.0);
    assert_eq!(spectrum.ms_level, 1);
    // A consumer that takes the source `PeakSpectrum` takes `MSSpectrum` here.
    let binned = BinnedSpectrum::new(&spectrum, BinConfig::default()).unwrap();
    assert!(binned.bins().is_empty());
    let boxed: Box<MSSpectrum> = Box::default();
    assert_eq!(*boxed, spectrum);
    drop(boxed);
}

/// `StandardTypes_test.cpp` `GOOD_TYPEDEF(PeakMap)` (declared twice in the
/// source): `PeakMap` is `MSExperiment`.
#[test]
fn peak_map_typedef_is_ms_experiment() {
    let experiment: MSExperiment = MSExperiment::new();
    assert_eq!(experiment, MSExperiment::default());
    assert!(experiment.spectra.is_empty());
    assert!(experiment.chromatograms.is_empty());
    let boxed: Box<MSExperiment> = Box::default();
    assert_eq!(*boxed, experiment);
    drop(boxed);
}

/// `Chromatogram` is `MSChromatogram`. The source test never instantiates this
/// alias (it repeats `PeakSpectrum` and `PeakMap` instead); this is an
/// additional native check of the third typedef.
#[test]
fn chromatogram_typedef_is_ms_chromatogram() {
    let chromatogram: MSChromatogram = MSChromatogram::new();
    assert_eq!(chromatogram, MSChromatogram::default());
    assert!(chromatogram.peaks.is_empty());
    assert_eq!(
        MSChromatogram::from_peaks(vec![ChromatogramPeak::new(1.0, 2.0)]).peaks[0],
        ChromatogramPeak {
            rt: 1.0,
            intensity: 2.0
        }
    );
    let boxed: Box<MSChromatogram> = Box::default();
    assert_eq!(*boxed, chromatogram);
    drop(boxed);
}
