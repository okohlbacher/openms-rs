// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! Section-by-section port of `BinnedSpectrum_test.cpp`
//! (source revision bc9cc12514c768385ce121d6ca4bb710fe1983c4, tier 3:
//! transcribed class-test literals, no C++ execution).
//!
//! The source fixture is `PILISSequenceDB_DFPIANGER_1.dta`, retained locally as
//! `tests/data/comparison_dfpianger.dta` (same bytes, see
//! `docs/COMPARISON_SUPPORT.md`). Its header line `1019.74 1` yields one
//! precursor at m/z 1019.74 with charge 1.

use openms::MSSpectrum;
use openms::comparison::{BinConfig, BinUnit, BinnedSpectrum};
use openms::format::dta;

/// `TEST_REAL_SIMILAR` default: absolute difference within 1e-5, or ratio
/// within 1 + 1e-5.
fn real_similar(a: f64, b: f64) {
    let absdiff = (a - b).abs();
    let ratio = if a.abs() > b.abs() { a / b } else { b / a };
    assert!(
        absdiff <= 1e-5 || (1.0 - ratio).abs() <= 1e-5,
        "{a} !~ {b} (absdiff {absdiff}, ratio {ratio})"
    );
}

fn s1() -> MSSpectrum {
    dta::read(include_bytes!("data/comparison_dfpianger.dta").as_slice()).unwrap()
}

/// Source `BinnedSpectrum(s1, 1.5, false, 2, 0.0)`.
fn absolute_1_5_spread_2() -> BinConfig {
    BinConfig {
        size: 1.5,
        unit: BinUnit::Absolute,
        spread: 2,
        offset: 0.0,
        ..BinConfig::default()
    }
}

/// Source `BinnedSpectrum(s1, 10, true, 0, 0.0)`: 10 ppm bins.
fn ppm_10() -> BinConfig {
    BinConfig {
        size: 10.0,
        unit: BinUnit::Ppm,
        spread: 0,
        offset: 0.0,
        ..BinConfig::default()
    }
}

/// Section `~BinnedSpectrum()`: the source deletes a null pointer; here a
/// constructed value is dropped.
#[test]
fn destructor() {
    let bs = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    drop(bs);
}

/// Section `BinnedSpectrum(const PeakSpectrum& ps, float size, UInt spread,
/// float offset)` (the section title omits the `bool unit_ppm` argument the
/// constructor actually takes).
#[test]
fn detailed_constructor() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    assert_eq!(bs1.config().size, 1.5);
    assert_eq!(bs1.config().spread, 2);
    assert_eq!(bs1.config().offset, 0.0);
    assert_eq!(bs1.config().unit, BinUnit::Absolute);
    assert!(!bs1.bins().is_empty());
}

/// Section `BinnedSpectrum(const BinnedSpectrum& source)`.
#[test]
fn copy_constructor() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    let copy = bs1.clone();
    assert_eq!(copy.config().size, bs1.config().size);
    assert_eq!(copy.precursors().len(), 1);
    assert_eq!(bs1.precursors().len(), 1);
    assert_eq!(
        copy.precursors()[0].mz as u32,
        bs1.precursors()[0].mz as u32
    );
    assert_eq!(copy.precursors()[0].mz as u32, 1019);
}

/// Section `BinnedSpectrum& operator=(const BinnedSpectrum& source)`: the
/// copy outlives a rebuilt original and still matches it.
#[test]
fn assignment_operator() {
    let mut bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    let copy = bs1.clone();
    bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    assert_eq!(copy.config().size, bs1.config().size);
    assert_eq!(
        copy.precursors()[0].mz as u32,
        bs1.precursors()[0].mz as u32
    );
    let mut assigned = BinnedSpectrum::new(&MSSpectrum::default(), ppm_10()).unwrap();
    assert!(assigned.bins().is_empty());
    assigned = copy.clone();
    assert_eq!(assigned, bs1);
}

/// Section `bool operator==(const BinnedSpectrum& rhs) const`.
#[test]
fn equality_operator() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    let copy = bs1.clone();
    assert!(bs1 == copy);
}

/// Section `bool operator!=(const BinnedSpectrum& rhs) const`.
#[test]
fn inequality_operator() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    let copy = bs1.clone();
    assert!(!(bs1 != copy));
}

/// Section `float getBinSize() const`.
#[test]
fn get_bin_size() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    assert_eq!(bs1.config().size, 1.5);
}

/// Section `UInt getBinSpread() const`.
#[test]
fn get_bin_spread() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    assert_eq!(bs1.config().spread, 2);
}

/// Section `SparseVectorIndexType getBinIndex(double mz) const` with 10 ppm
/// bins.
#[test]
fn get_bin_index() {
    let bs1 = BinnedSpectrum::new(&s1(), ppm_10()).unwrap();
    let cfg = bs1.config();
    assert_eq!(cfg.bin_index(1.0).unwrap(), 0);
    assert_eq!(cfg.bin_index(10.0).unwrap(), 230259);
    assert_eq!(cfg.bin_index(100.0).unwrap(), 460519);
    assert_eq!(cfg.bin_index(1000.0).unwrap(), 690778);
}

/// Section `float getBinLowerMZ(size_t i) const`: all nineteen source
/// assertions, ppm bins first, then 1.0 m/z bins with 0.5 offset.
#[test]
fn get_bin_lower_mz() {
    let bs1 = BinnedSpectrum::new(&s1(), ppm_10()).unwrap();
    let cfg = bs1.config();
    let lower = |i: usize| f64::from(cfg.bin_lower_mz(i).unwrap());
    let index = |mz: f64| cfg.bin_index(mz).unwrap();
    real_similar(lower(0), 1.0); // m/z = 1 corresponds to lowest index
    real_similar(lower(1), 1.0 + 10.0 * 1e-6); // (1 + 10 ppm)
    real_similar(lower(1000), (1.0_f64 + 10.0 * 1e-6).powi(1000));
    real_similar(lower(index(1.0)), 1.0);
    real_similar(lower(index(10.0)), 10.0);
    real_similar(lower(index(100.0)), 100.0);
    real_similar(lower(index(1000.0)), 1000.0);
    drop(bs1);

    // 1.0 m/z bins with 0.5 offset: floats close to nominal masses fall into
    // the same bin.
    let bs2 = BinnedSpectrum::new(
        &s1(),
        BinConfig {
            size: 1.0,
            unit: BinUnit::Absolute,
            spread: 0,
            offset: 0.5,
            ..BinConfig::default()
        },
    )
    .unwrap();
    let cfg = bs2.config();
    let lower = |i: usize| f64::from(cfg.bin_lower_mz(i).unwrap());
    let index = |mz: f64| cfg.bin_index(mz).unwrap();
    assert_eq!(index(999.99), index(1000.01));
    assert_eq!(index(99.99), index(100.01));
    assert_eq!(index(9.99), index(10.01));
    assert_eq!(index(0.99), index(1.01));
    assert_eq!(index(0.0), index(0.01));
    real_similar(lower(0), -0.5); // because of offset, bin starts at -0.5
    real_similar(lower(1), 0.5);
    real_similar(lower(1000), 999.5);
    real_similar(lower(index(0.5)), 0.5);
    real_similar(lower(index(9.5)), 9.5);
    real_similar(lower(index(99.5)), 99.5);
    real_similar(lower(index(999.5)), 999.5);
}

/// Section `const SparseVectorType* getBins() const`: 347 stored bins, bin
/// 658 holds 501645, and a read never changes the stored count.
#[test]
fn get_bins_const() {
    let bs1 = BinnedSpectrum::new(
        &s1(),
        BinConfig {
            offset: 0.0,
            ..absolute_1_5_spread_2()
        },
    )
    .unwrap();
    assert_eq!(bs1.bins().len(), 347);
    assert_eq!(bs1.bins()[&658], 501645.0);
    assert_eq!(bs1.bins().len(), 347);
    let c = bs1.bins().iter().count();
    assert_eq!(bs1.bins().len(), c);
    // Source coeffRef on a missing index inserts a zero; the native lookup
    // does not, so the count stays 347 even after a miss.
    assert_eq!(bs1.bins().get(&0), None);
    assert_eq!(bs1.bins().len(), 347);
}

/// Section `SparseVectorType* getBins()`: the source reads through the
/// mutable accessor; the native container is read-only and the m/z lookup
/// returns the same coefficient. Bin 658 with size 1.5 covers [987, 988.5).
#[test]
fn get_bins_mutable_equivalent() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    assert_eq!(bs1.bins()[&658], 501645.0);
    assert_eq!(bs1.bin_intensity(987.0).unwrap(), 501645.0);
    assert_eq!(bs1.config().bin_index(987.0).unwrap(), 658);
    real_similar(f64::from(bs1.config().bin_lower_mz(658).unwrap()), 987.0);
}

/// Section `void setBinning()`: `NOT_TESTABLE` in the source (private
/// `binSpectrum_`, exercised by construction). Native check: spread 2 adds
/// each peak to two neighbours on either side, stopping at bin 0.
#[test]
fn set_binning_is_construction() {
    let spectrum = MSSpectrum::from_peaks(vec![openms::Peak1D::new(0.0, 4.0)]);
    let bs = BinnedSpectrum::new(&spectrum, absolute_1_5_spread_2()).unwrap();
    assert_eq!(
        bs.bins().iter().map(|(&i, &v)| (i, v)).collect::<Vec<_>>(),
        [(0, 4.0), (1, 4.0), (2, 4.0)]
    );
}

/// Section `bool BinnedSpectrum::isCompatible(const BinnedSpectrum& a, const
/// BinnedSpectrum& b)`.
#[test]
fn is_compatible() {
    let bs1 = BinnedSpectrum::new(&s1(), absolute_1_5_spread_2()).unwrap();
    let bs2 = BinnedSpectrum::new(
        &s1(),
        BinConfig {
            size: 1.234,
            ..absolute_1_5_spread_2()
        },
    )
    .unwrap();
    assert!(!bs1.is_compatible(&bs2));
    assert!(bs1.is_compatible(&bs1));
}

/// Source class constants and their relation to the `BinConfig` defaults.
#[test]
fn recommended_layout_constants() {
    assert_eq!(BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES, 1.0005);
    assert_eq!(BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES, 0.02);
    assert_eq!(BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES, 0.0);
    assert_eq!(BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES, 0.4);
    let default = BinConfig::default();
    assert_eq!(default.size, BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES);
    assert_eq!(default.offset, BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES);
    assert_eq!(default.spread, 0);
}
