// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
// $Authors: Chris Bielow, OpenMS Rust contributors $

//! The signal-to-noise estimator base and the random-scan noise estimate.
//!
//! Ports `src/openms/include/OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h`
//! and `src/openms/source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp`.
//! `docs/SIGNAL_TO_NOISE_SUPPORT.md` holds the API mapping, the preserved source
//! conventions, the native differences and the evidence.
//!
//! The source header is the abstract base class of a signal-to-noise
//! estimator: an estimator provides the signal-to-noise ratio of every raw data
//! point of a container. Its two concrete subclasses in the pinned core are
//! ported as
//! [`SignalToNoiseEstimatorMedian`](crate::processing::peak_picking::SignalToNoiseEstimatorMedian)
//! and
//! [`SignalToNoiseEstimatorMeanIterative`](crate::processing::mean_noise::SignalToNoiseEstimatorMeanIterative),
//! and both implement [`SignalToNoiseEstimator`], the counterpart of the pure
//! virtual `computeSTN_` contract.
//!
//! # Linux x86-64 Release arithmetic
//!
//! The port follows the Linux x86-64 Release build of the pinned core wherever
//! the platform decides a result. The functions of the private `x86` module
//! reproduce, bit for bit, the instructions that build emits at the sites the
//! support document lists: the SSE rules for NaN operands (the first operand's
//! NaN wins, an invalid operation yields the negative default NaN) and the
//! truncating conversions `cvttsd2si`, whose out-of-range and NaN result is the
//! "integer indefinite" value. Each use names the source line and the
//! instruction it emulates.

use crate::kernel::{MSChromatogram, MSExperiment, MSSpectrum};
use crate::processing::mean_noise::{MeanNoiseEstimates, SignalToNoiseEstimatorMeanIterative};
use crate::{Error, Result};

/// The contract of a signal-to-noise estimator: the ratio of every data point
/// of a container.
///
/// Source abstract class `SignalToNoiseEstimator<Container>`. Its pure virtual
/// `computeSTN_(const Container&)`, called by `init(container)`, fills the
/// member `stn_estimates_`, which `getSignalToNoise(index)` reads. Here an
/// estimator is stateless: [`compute_stn`](Self::compute_stn) returns the
/// complete estimation as a value, and
/// [`signal_to_noise`](Self::signal_to_noise) is the vector the source serves
/// one index at a time.
///
/// The source template parameter `Container` (default `MSSpectrum`) and its
/// typedefs `PeakIterator` and `PeakType` become the separate entry points for
/// slices, spectra and chromatograms.
///
/// The ported implementations use their native input contract in these trait
/// methods; each type documents the inherent method that selects the source
/// behaviour.
pub trait SignalToNoiseEstimator {
    /// One estimation: the ratios plus whatever the estimator reports besides.
    type Estimates;

    /// Source `init(container)`, i.e. `computeSTN_`, over parallel position and
    /// intensity slices.
    ///
    /// # Errors
    ///
    /// The implementation's input and option refusals.
    fn compute_stn(&self, positions: &[f64], intensities: &[f64]) -> Result<Self::Estimates>;

    /// Source `getSignalToNoise(index)` for every index of an estimation, in
    /// input order.
    ///
    /// The source checks the index only with `OPENMS_POSTCONDITION`, which a
    /// Release build compiles out, and then reads past the end of its vector;
    /// a slice's `get` returns `None` there instead.
    fn signal_to_noise(estimates: &Self::Estimates) -> &[f64];

    /// Source `init(spectrum)` for `Container = MSSpectrum`: positions are m/z,
    /// intensities the stored `f32` values.
    ///
    /// # Errors
    ///
    /// As [`compute_stn`](Self::compute_stn).
    fn compute_stn_spectrum(&self, spectrum: &MSSpectrum) -> Result<Self::Estimates> {
        let positions: Vec<f64> = spectrum.peaks.iter().map(|p| p.mz).collect();
        let intensities: Vec<f64> = spectrum
            .peaks
            .iter()
            .map(|p| x86::widen(p.intensity))
            .collect();
        self.compute_stn(&positions, &intensities)
    }

    /// Source `init(chromatogram)` for `Container = MSChromatogram`: positions
    /// are retention times in seconds.
    ///
    /// # Errors
    ///
    /// As [`compute_stn`](Self::compute_stn).
    fn compute_stn_chromatogram(&self, chromatogram: &MSChromatogram) -> Result<Self::Estimates> {
        let positions: Vec<f64> = chromatogram.peaks.iter().map(|p| p.rt).collect();
        let intensities: Vec<f64> = chromatogram
            .peaks
            .iter()
            .map(|p| x86::widen(p.intensity))
            .collect();
        self.compute_stn(&positions, &intensities)
    }
}

impl SignalToNoiseEstimator for SignalToNoiseEstimatorMeanIterative {
    type Estimates = MeanNoiseEstimates;

    /// [`SignalToNoiseEstimatorMeanIterative::estimate`].
    fn compute_stn(&self, positions: &[f64], intensities: &[f64]) -> Result<MeanNoiseEstimates> {
        self.estimate(positions, intensities)
    }

    fn signal_to_noise(estimates: &MeanNoiseEstimates) -> &[f64] {
        &estimates.signal_to_noise
    }

    /// [`SignalToNoiseEstimatorMeanIterative::estimate_spectrum`], which also
    /// validates the spectrum.
    fn compute_stn_spectrum(&self, spectrum: &MSSpectrum) -> Result<MeanNoiseEstimates> {
        self.estimate_spectrum(spectrum)
    }

    /// [`SignalToNoiseEstimatorMeanIterative::estimate_chromatogram`], which
    /// also validates the chromatogram.
    fn compute_stn_chromatogram(
        &self,
        chromatogram: &MSChromatogram,
    ) -> Result<MeanNoiseEstimates> {
        self.estimate_chromatogram(chromatogram)
    }
}

/// Mean and variance of a Gaussian fitted to intensities.
///
/// Source protected struct `SignalToNoiseEstimator::GaussianEstimate`
/// ("parameters my, sigma for a Gaussian distribution", accessors `mean` and
/// `variance`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussianEstimate {
    /// Mean of the estimated Gaussian.
    pub mean: f64,
    /// Variance of the estimated Gaussian (the population variance, divided by
    /// `n`).
    pub variance: f64,
}

impl GaussianEstimate {
    /// Mean and population variance of `values`, as source protected
    /// `estimate_(first, last)` computes them: one pass summing in input order,
    /// a division by the count, a second pass summing the squared deviations
    /// `(mean - value)^2` in input order and a division by the count.
    ///
    /// An empty slice divides zero by zero, so both fields are the Linux x86-64
    /// default NaN (bits `0xfff8000000000000`), as in the source. NaN operands
    /// propagate by the SSE rules the Release build follows (see the module
    /// documentation).
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::processing::noise_estimation::GaussianEstimate;
    ///
    /// let g = GaussianEstimate::of(&[1.0, 2.0, 3.0, 6.0]);
    /// assert_eq!(g.mean, 3.0);
    /// assert_eq!(g.variance, 3.5);
    /// assert_eq!(GaussianEstimate::of(&[]).mean.to_bits(), 0xfff8_0000_0000_0000);
    /// ```
    pub fn of(values: &[f64]) -> Self {
        Self::of_indexed(values.len(), |i| values[i])
    }

    /// [`GaussianEstimate::of`] over `n` values read by index.
    pub(crate) fn of_indexed(n: usize, value: impl Fn(usize) -> f64) -> Self {
        // SignalToNoiseEstimator.h:115-137. The Release build keeps the running
        // sum as the first operand (`addsd I, m`, `subsd I, tmp`, `addsd sq, v`)
        // and converts the count with `cvtsi2sd` from a 32-bit int.
        let mut mean = 0.0;
        for i in 0..n {
            mean = x86::add(mean, value(i));
        }
        mean = x86::div(mean, n as f64);
        let mut variance = 0.0;
        for i in 0..n {
            let deviation = x86::sub(mean, value(i));
            variance = x86::add(variance, x86::mul(deviation, deviation));
        }
        variance = x86::div(variance, n as f64);
        Self { mean, variance }
    }
}

/// The source's `std::default_random_engine`: libstdc++'s `minstd_rand0`, with
/// the draws of `std::uniform_real_distribution<double>(0, 1)` as the Release
/// build computes them.
///
/// A Lehmer generator with multiplier `16807` and modulus `2^31 - 1`
/// (`linear_congruential_engine<uint_fast32_t, 16807, 0, 2147483647>`, the
/// typedef at `bits/random.h:1743` of GCC 14.4). Seeding keeps
/// `seed mod (2^31 - 1)` and replaces a zero state by one
/// (`bits/random.tcc:119-126`).
///
/// # Examples
///
/// ```
/// use openms::processing::noise_estimation::MinstdRand0;
///
/// let mut engine = MinstdRand0::new(1);
/// assert_eq!(engine.next_u32(), 16807);
/// assert_eq!(engine.next_u32(), 282_475_249);
/// ```
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MinstdRand0 {
    state: u64,
}

impl MinstdRand0 {
    /// The multiplier, `16807`.
    pub const MULTIPLIER: u64 = 16_807;
    /// The modulus, `2^31 - 1`.
    pub const MODULUS: u64 = 2_147_483_647;

    /// An engine seeded with `seed`, the value the source passes after
    /// converting `time(nullptr)` to the engine's `unsigned long` result type
    /// (a negative time wraps modulo `2^64`).
    pub fn new(seed: u64) -> Self {
        let state = seed % Self::MODULUS;
        Self {
            state: if state == 0 { 1 } else { state },
        }
    }

    /// The next engine output, in `1..2^31 - 1`.
    pub fn next_u32(&mut self) -> u32 {
        // state < 2^31 and the multiplier < 2^15, so the product fits u64.
        self.state = self.state * Self::MULTIPLIER % Self::MODULUS;
        u32::try_from(self.state).unwrap_or(u32::MAX)
    }

    /// A draw of `std::uniform_real_distribution<double>(0.0, 1.0)`, in
    /// `[0, 1)`.
    ///
    /// `std::generate_canonical<double, 53>` (`bits/random.tcc:3349-3381`)
    /// consumes two engine outputs, `a` then `b`, and returns
    /// `((b - 1) * r + (a - 1)) / r^2` with `r = 2^31 - 2`, clamped below one
    /// to `nextafter(1, 0)`. The Release build folds `r^2` to the double
    /// `4611686009837453312` (`0x43cfffffff000000`), which is `r^2` rounded, and
    /// adds `a - 1` last (`SignalToNoiseEstimator.cpp:42` in
    /// `libOpenMS.so` at `0x18684f6`).
    pub fn uniform01(&mut self) -> f64 {
        const R: f64 = 2_147_483_646.0;
        const R_SQUARED: f64 = 4_611_686_009_837_453_312.0;
        const BELOW_ONE: f64 = 0.999_999_999_999_999_9;
        let first = f64::from(self.next_u32() - 1);
        let second = f64::from(self.next_u32() - 1);
        let sum = second * R + (first + 0.0);
        let value = sum / R_SQUARED;
        if value >= 1.0 { BELOW_ONE } else { value }
    }
}

/// Native ceiling on the intensities [`estimate_noise_from_random_scans`]
/// copies and selects over, across all drawn scans.
pub const RANDOM_SCAN_NOISE_MAX_WORK: usize = 100_000_000;

/// Options of the random-scan noise estimate.
///
/// Source free function `estimateNoiseFromRandomScans(exp, ms_level, n_scans =
/// 10, percentile = 80)`, whose random engine is seeded with `time(nullptr)`.
/// The seed is explicit here, which the crate's determinism contract requires.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RandomScanNoise {
    /// MS level of the spectra that count as candidates.
    pub ms_level: u32,
    /// Number of scans to draw, source default `10`.
    pub n_scans: u32,
    /// Percentile of each drawn scan's intensities, source default `80`.
    pub percentile: f64,
    /// Engine seed; the source uses `time(nullptr)` converted to `unsigned
    /// long`.
    pub seed: u64,
    /// Native ceiling on copied and selected intensities plus drawn scans,
    /// default [`RANDOM_SCAN_NOISE_MAX_WORK`]. The source has none.
    pub max_work: usize,
}

impl RandomScanNoise {
    /// The source defaults for `ms_level` with an explicit `seed`.
    pub fn new(ms_level: u32, seed: u64) -> Self {
        Self {
            ms_level,
            n_scans: 10,
            percentile: 80.0,
            seed,
            max_work: RANDOM_SCAN_NOISE_MAX_WORK,
        }
    }

    /// The average intensity at `percentile` of `n_scans` randomly drawn scans.
    ///
    /// Source `estimateNoiseFromRandomScans` (`SignalToNoiseEstimator.cpp:21-53`),
    /// reproduced with its defined quirks:
    ///
    /// * The candidates are the non-empty spectra of `ms_level`, but the drawn
    ///   index addresses the **experiment**, not the candidate list (`:44`
    ///   reads `exp[scan]`), so the level filter is discarded and a drawn
    ///   spectrum may have any level.
    /// * The index is `(UInt)(u * (candidates - 1))` (`:42`): the last candidate
    ///   position is never drawn, and one candidate always draws experiment
    ///   index `0`. The Release build truncates with a 64-bit `cvttsd2si` and
    ///   keeps the low 32 bits.
    /// * The position in a drawn scan is `(Size)(size * percentile / 100.0)`
    ///   (`:48`), converted as the Release build converts a double to
    ///   `unsigned long`; a value in `(-1, 0)` selects the minimum, and a value
    ///   of `2^64` or more, including infinity, selects index `0` (the
    ///   `subsd 2^63; cvttsd2si; btc 63` path).
    /// * The selected intensity is the one `std::nth_element` puts at that
    ///   position, summed in `f32` in draw order and divided by `n_scans` in
    ///   `f32` (`:50-52`). Without NaN intensities that is the order statistic
    ///   (a `-0.0` and a `+0.0` are interchangeable there, and the sum, which
    ///   starts at `+0`, cannot tell them apart). With NaN intensities
    ///   `operator<` is not a strict weak ordering, which the standard makes
    ///   undefined; the Release toolchain's `nth_element` nevertheless stays in
    ///   bounds and is deterministic for every irreflexive, asymmetric
    ///   comparison, which `<` on `float` is, so this port runs that exact
    ///   algorithm (see the private `libstdcxx` module) and returns what that
    ///   build returns.
    /// * `n_scans = 0` divides zero by zero and returns the negative default
    ///   `f32` NaN.
    /// * Without candidates the result is `0.0` (`:32`); the engine is not
    ///   seeded then, as the source does not call `time` there.
    ///
    /// # Errors
    ///
    /// [`Error::Unsupported`] exactly where the source becomes undefined, at
    /// the draw where it happens:
    ///
    /// * the position exceeds the drawn scan's size, so `tmp.begin() + idx` at
    ///   `:49` points past the end (a percentile above `100`, a percentile of
    ///   `-100 / size` or below, or a NaN product);
    /// * the position equals the drawn scan's size, so `tmp[idx]` at `:50`
    ///   reads past the end (a percentile of exactly `100`, or an empty drawn
    ///   scan with a percentile whose product with zero is not NaN).
    ///
    /// [`Error::InvalidValue`] when the native `max_work` ceiling would be
    /// exceeded.
    ///
    /// # Examples
    ///
    /// ```
    /// use openms::kernel::{MSExperiment, MSSpectrum, Peak1D};
    /// use openms::processing::noise_estimation::RandomScanNoise;
    ///
    /// let mut experiment = MSExperiment::new();
    /// for level in [2, 1] {
    ///     let mut spectrum = MSSpectrum::from_peaks(
    ///         (0..5).map(|i| Peak1D::new(f64::from(i), (i * i) as f32)).collect(),
    ///     );
    ///     spectrum.ms_level = level;
    ///     experiment.spectra.push(spectrum);
    /// }
    /// // One MS1 candidate: every draw reads experiment index 0, an MS2 scan.
    /// let noise = RandomScanNoise::new(1, 42).estimate(&experiment).unwrap();
    /// assert_eq!(noise, 16.0);
    /// ```
    pub fn estimate(&self, experiment: &MSExperiment) -> Result<f32> {
        let candidates = experiment
            .spectra
            .iter()
            .filter(|s| s.ms_level == self.ms_level && !s.peaks.is_empty())
            .count();
        if candidates == 0 {
            return Ok(0.0);
        }
        let mut engine = MinstdRand0::new(self.seed);
        if self.n_scans == 0 {
            // `noise / n_scans` with both zero: `divss` of 0 by 0.
            return Ok(x86::div32(0.0, 0.0));
        }
        let scale = (candidates - 1) as f64;
        let mut noise = 0.0_f32;
        let mut work = 0usize;
        let mut values: Vec<f32> = Vec::new();
        for draw in 0..self.n_scans {
            let u = engine.uniform01();
            // `(UInt)(u * (size - 1))`: 64-bit truncation, low 32 bits kept.
            let scan = x86::cvttsd2si64((u + 0.0) * scale) as u32;
            let spectrum = experiment.spectra.get(scan as usize).ok_or_else(|| {
                Error::InvalidValue(format!(
                    "random-scan index {scan} is outside the experiment"
                ))
            })?;
            work = work
                .checked_add(spectrum.peaks.len())
                .and_then(|w| w.checked_add(1))
                .filter(|w| *w <= self.max_work)
                .ok_or_else(|| {
                    Error::InvalidValue(
                        "random-scan noise estimation exceeds configured work limit".into(),
                    )
                })?;
            values.clear();
            values.extend(spectrum.peaks.iter().map(|p| p.intensity));
            let len = values.len();
            let index = x86::f64_to_u64(x86::div(x86::mul(len as f64, self.percentile), 100.0));
            let context = |what: &str| {
                Error::Unsupported(format!(
                    "estimateNoiseFromRandomScans is undefined here: draw {draw} reads experiment \
                     index {scan} with {len} intensities and position {index}; {what}"
                ))
            };
            if index > len as u64 {
                return Err(context(
                    "SignalToNoiseEstimator.cpp:49 advances the iterator past the end",
                ));
            }
            if index == len as u64 {
                return Err(context(
                    "SignalToNoiseEstimator.cpp:50 reads one element past the end",
                ));
            }
            // index < len, so it fits usize.
            let position = usize::try_from(index).unwrap_or(0);
            // :49, with the Release toolchain's own algorithm, so that NaN
            // intensities land where that build puts them.
            libstdcxx::nth_element(&mut values, position);
            // `addss noise, tmp[idx]`: the selected value is the first operand.
            noise = x86::add32(values[position], noise);
        }
        Ok(x86::div32(noise, self.n_scans as f32))
    }
}

/// Source free function `estimateNoiseFromRandomScans(exp, ms_level, n_scans,
/// percentile)` with an explicit `seed` in place of `time(nullptr)`: "picks
/// `n_scans` from the given `ms_level` randomly and returns either average
/// intensity at a certain `percentile`. If no scans with the required level
/// are present, 0.0 is returned". See [`RandomScanNoise::estimate`] for the
/// exact behaviour, which differs from that description.
///
/// # Errors
///
/// As [`RandomScanNoise::estimate`], with the default work ceiling.
pub fn estimate_noise_from_random_scans(
    experiment: &MSExperiment,
    ms_level: u32,
    n_scans: u32,
    percentile: f64,
    seed: u64,
) -> Result<f32> {
    RandomScanNoise {
        n_scans,
        percentile,
        ..RandomScanNoise::new(ms_level, seed)
    }
    .estimate(experiment)
}

/// `std::nth_element` over `float` values with `operator<`, as the Release
/// toolchain's C++ library implements it.
///
/// A line-by-line port of `std::nth_element`, `__introselect`,
/// `__unguarded_partition_pivot`, `__move_median_to_first`,
/// `__unguarded_partition`, `__insertion_sort`, `__unguarded_linear_insert`
/// and `__heap_select` from `bits/stl_algo.h` (sha256 `0598c5b1…`), and
/// `__make_heap`, `__adjust_heap`, `__pop_heap` and `__push_heap` from
/// `bits/stl_heap.h` (sha256 `f18f83b2…`) of the conda-forge GCC 14.4.0 that
/// built `libOpenMS.so`; `__lg` is `bit_width(n) - 1` (`bits/stl_algobase.h`).
/// Every comparison is made in the library's order and with its operands, and
/// every move is the library's, so the whole permutation, not only the
/// selected value, is the one that build produces.
///
/// # Memory safety for every `float` input
///
/// The source calls this with NaN intensities too, where `<` is not a strict
/// weak ordering. The algorithm still stays in bounds, because its
/// "unguarded" loops rely only on comparisons it has made, plus the
/// irreflexivity (`!(x < x)`) and asymmetry (`x < y` implies `!(y < x)`) that
/// `<` on `float` keeps with NaN:
///
/// * `__move_median_to_first(first, a, b, c)` leaves, in each of its six
///   branches, one of the two samples it did not move that is not less than
///   the pivot, established by one of its comparisons (by asymmetry where
///   that comparison found the pivot less than it).
/// * The upward scan of `__unguarded_partition` therefore stops at that
///   sample on its first pass; on later passes it stops at the element the
///   previous swap put at `last`, which the scan had found not less than the
///   pivot. The downward scan stops at the pivot's own slot at the latest,
///   because the pivot is not less than itself, and the pivot does not move
///   during the partition. The returned cut lies in `first + 1 .. last`, so
///   `__introselect`'s range shrinks on every pass.
/// * `__unguarded_linear_insert` is called for an element found not less
///   than the range's first element, which it does not move, so the scan
///   stops there at the latest.
/// * The heap functions index by explicit lengths only.
///
/// The Rust loops below carry explicit bounds that, by this argument, never
/// decide anything.
pub(crate) mod libstdcxx {
    /// `__gnu_cxx::__ops::__iter_less_iter` and friends: `lhs < rhs`.
    fn lt(a: f32, b: f32) -> bool {
        a < b
    }

    /// `std::nth_element(v.begin(), v.begin() + nth, v.end())`. An `nth` of
    /// `v.len()` or more returns at once, as the library does for
    /// `nth == last` (callers never pass more).
    pub(crate) fn nth_element(v: &mut [f32], nth: usize) {
        let n = v.len();
        if n == 0 || nth >= n {
            return;
        }
        // `std::__lg(last - first) * 2`; n > 0.
        let lg = (usize::BITS - 1 - n.leading_zeros()) as usize;
        introselect(v, 0, nth, n, lg * 2);
    }

    fn introselect(v: &mut [f32], mut first: usize, nth: usize, mut last: usize, mut depth: usize) {
        while last - first > 3 {
            if depth == 0 {
                heap_select(v, first, nth + 1, last);
                // Place the nth largest element in its final position.
                v.swap(first, nth);
                return;
            }
            depth -= 1;
            let cut = unguarded_partition_pivot(v, first, last);
            if cut <= nth {
                first = cut;
            } else {
                last = cut;
            }
        }
        insertion_sort(v, first, last);
    }

    fn unguarded_partition_pivot(v: &mut [f32], first: usize, last: usize) -> usize {
        let mid = first + (last - first) / 2;
        move_median_to_first(v, first, first + 1, mid, last - 1);
        unguarded_partition(v, first + 1, last, first)
    }

    /// Swaps the median of `v[a]`, `v[b]` and `v[c]` under `<` into `result`.
    fn move_median_to_first(v: &mut [f32], result: usize, a: usize, b: usize, c: usize) {
        let median = if lt(v[a], v[b]) {
            if lt(v[b], v[c]) {
                b
            } else if lt(v[a], v[c]) {
                c
            } else {
                a
            }
        } else if lt(v[a], v[c]) {
            a
        } else if lt(v[b], v[c]) {
            c
        } else {
            b
        };
        v.swap(result, median);
    }

    fn unguarded_partition(
        v: &mut [f32],
        mut first: usize,
        mut last: usize,
        pivot: usize,
    ) -> usize {
        let end = last;
        loop {
            while first < end && lt(v[first], v[pivot]) {
                first += 1;
            }
            last -= 1;
            while last > pivot && lt(v[pivot], v[last]) {
                last -= 1;
            }
            if first >= last {
                return first;
            }
            v.swap(first, last);
            first += 1;
        }
    }

    fn insertion_sort(v: &mut [f32], first: usize, last: usize) {
        if first == last {
            return;
        }
        for i in first + 1..last {
            if lt(v[i], v[first]) {
                let value = v[i];
                // _GLIBCXX_MOVE_BACKWARD3(first, i, i + 1)
                v.copy_within(first..i, first + 1);
                v[first] = value;
            } else {
                unguarded_linear_insert(v, first, i);
            }
        }
    }

    /// `__unguarded_linear_insert(last)` with `__val_comp_iter`; `floor` is
    /// the range's first element, which the caller found not greater.
    fn unguarded_linear_insert(v: &mut [f32], floor: usize, mut last: usize) {
        let value = v[last];
        let mut next = last - 1;
        while lt(value, v[next]) {
            v[last] = v[next];
            last = next;
            if next == floor {
                break;
            }
            next -= 1;
        }
        v[last] = value;
    }

    fn heap_select(v: &mut [f32], first: usize, middle: usize, last: usize) {
        make_heap(v, first, middle);
        for i in middle..last {
            if lt(v[i], v[first]) {
                pop_heap(v, first, middle, i);
            }
        }
    }

    fn make_heap(v: &mut [f32], first: usize, last: usize) {
        if last - first < 2 {
            return;
        }
        let len = last - first;
        let mut parent = (len - 2) / 2;
        loop {
            let value = v[first + parent];
            adjust_heap(v, first, parent, len, value);
            if parent == 0 {
                return;
            }
            parent -= 1;
        }
    }

    fn adjust_heap(v: &mut [f32], first: usize, mut hole: usize, len: usize, value: f32) {
        let top = hole;
        let mut second = hole;
        while second < (len - 1) / 2 {
            second = 2 * (second + 1);
            if lt(v[first + second], v[first + second - 1]) {
                second -= 1;
            }
            v[first + hole] = v[first + second];
            hole = second;
        }
        if len & 1 == 0 && second == (len - 2) / 2 {
            second = 2 * (second + 1);
            v[first + hole] = v[first + second - 1];
            hole = second - 1;
        }
        push_heap(v, first, hole, top, value);
    }

    fn pop_heap(v: &mut [f32], first: usize, last: usize, result: usize) {
        let value = v[result];
        v[result] = v[first];
        adjust_heap(v, first, 0, last - first, value);
    }

    /// `__push_heap` with `__iter_comp_val`; `(hole - 1) / 2` truncates
    /// toward zero in the signed source, so a hole of `0` has parent `0`.
    fn push_heap(v: &mut [f32], first: usize, mut hole: usize, top: usize, value: f32) {
        let mut parent = hole.saturating_sub(1) / 2;
        while hole > top && lt(v[first + parent], value) {
            v[first + hole] = v[first + parent];
            hole = parent;
            parent = hole.saturating_sub(1) / 2;
        }
        v[first + hole] = value;
    }

    #[cfg(test)]
    mod tests {
        use super::nth_element;
        use std::collections::BTreeMap;

        /// Tier 1: the whole vector `std::nth_element` leaves, printed by the
        /// Release toolchain for the 390 vectors of `nth_vectors.tsv`
        /// (oracle case `nth_vectors`, `../oracle/sne-completion/`), 30 of
        /// which reach the `__heap_select` fallback, compared bit for bit.
        #[test]
        fn permutations_match_the_release_toolchain() {
            let oracle: BTreeMap<&str, Vec<u32>> =
                include_str!("../../tests/data/signal_to_noise/oracle.tsv")
                    .lines()
                    .filter(|l| l.starts_with("nthout\t"))
                    .map(|l| {
                        let c: Vec<&str> = l.split('\t').collect();
                        let bits = c[3..]
                            .iter()
                            .map(|h| u32::from_str_radix(h, 16).unwrap())
                            .collect();
                        (c[2], bits)
                    })
                    .collect();
            let mut compared = 0;
            for line in include_str!("../../tests/data/signal_to_noise/nth_vectors.tsv").lines() {
                let c: Vec<&str> = line.split('\t').collect();
                let nth: usize = c[1].parse().unwrap();
                let mut values: Vec<f32> = c[2]
                    .split(',')
                    .map(|h| f32::from_bits(u32::from_str_radix(h, 16).unwrap()))
                    .collect();
                nth_element(&mut values, nth);
                let bits: Vec<u32> = values.iter().map(|v| v.to_bits()).collect();
                assert_eq!(Some(&bits), oracle.get(c[0]), "{}", c[0]);
                compared += 1;
            }
            assert_eq!((compared, oracle.len()), (390, 390));
        }

        #[test]
        fn selects_the_order_statistic_without_nan() {
            // Independent check: without NaN the selected value is the sorted
            // value, for every position and for sizes on both sides of the
            // insertion-sort threshold and past the depth limit.
            let mut state = 12_345_u32;
            for n in 1..80usize {
                for trial in 0..8 {
                    let values: Vec<f32> = (0..n)
                        .map(|_| {
                            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                            ((state >> 16) % (3 + trial as u32 * 7)) as f32
                        })
                        .collect();
                    let mut sorted = values.clone();
                    sorted.sort_by(f32::total_cmp);
                    for nth in 0..n {
                        let mut v = values.clone();
                        nth_element(&mut v, nth);
                        assert_eq!(v[nth], sorted[nth], "n {n} trial {trial} nth {nth}");
                        assert!(v[..nth].iter().all(|x| *x <= v[nth]));
                        assert!(v[nth..].iter().all(|x| *x >= v[nth]));
                        let mut back = v.clone();
                        back.sort_by(f32::total_cmp);
                        assert_eq!(back, sorted);
                    }
                }
            }
        }

        #[test]
        fn nan_inputs_stay_in_bounds_and_keep_the_multiset() {
            let pool = [
                f32::NAN,
                -f32::NAN,
                0.0,
                -0.0,
                1.0,
                2.0,
                f32::INFINITY,
                -1.0,
            ];
            let mut state = 7_u32;
            for n in 1..70usize {
                for _ in 0..20 {
                    let values: Vec<f32> = (0..n)
                        .map(|_| {
                            state = state.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                            pool[(state >> 16) as usize % pool.len()]
                        })
                        .collect();
                    for nth in [0, n / 2, n - 1] {
                        let mut v = values.clone();
                        nth_element(&mut v, nth);
                        let mut a: Vec<u32> = v.iter().map(|x| x.to_bits()).collect();
                        let mut b: Vec<u32> = values.iter().map(|x| x.to_bits()).collect();
                        a.sort_unstable();
                        b.sort_unstable();
                        assert_eq!(a, b);
                    }
                }
            }
        }
    }
}

/// Bit-exact emulations of the Linux x86-64 Release build's scalar floating
/// point instructions where the platform decides the result.
pub(crate) mod x86 {
    /// The SSE "QNaN floating-point indefinite": what an invalid operation on
    /// non-NaN operands returns.
    pub(crate) const DEFAULT_NAN: u64 = 0xfff8_0000_0000_0000;
    /// The single-precision indefinite.
    pub(crate) const DEFAULT_NAN_F32: u32 = 0xffc0_0000;
    const QUIET: u64 = 0x0008_0000_0000_0000;
    const QUIET_F32: u32 = 0x0040_0000;

    fn default_nan() -> f64 {
        f64::from_bits(DEFAULT_NAN)
    }

    /// A binary SSE operation (`addsd`, `subsd`, `mulsd`, `divsd`) whose first
    /// source operand is `a`: a NaN operand is returned quieted, `a`'s first;
    /// an invalid operation returns the default NaN.
    fn binary(a: f64, b: f64, result: f64) -> f64 {
        if a.is_nan() {
            f64::from_bits(a.to_bits() | QUIET)
        } else if b.is_nan() {
            f64::from_bits(b.to_bits() | QUIET)
        } else if result.is_nan() {
            default_nan()
        } else {
            result
        }
    }

    /// `addsd`: `a + b` with `a` the first operand.
    pub(crate) fn add(a: f64, b: f64) -> f64 {
        binary(a, b, a + b)
    }

    /// `subsd`: `a - b`.
    pub(crate) fn sub(a: f64, b: f64) -> f64 {
        binary(a, b, a - b)
    }

    /// `mulsd`: `a * b` with `a` the first operand.
    pub(crate) fn mul(a: f64, b: f64) -> f64 {
        binary(a, b, a * b)
    }

    /// `divsd`: `a / b`.
    pub(crate) fn div(a: f64, b: f64) -> f64 {
        binary(a, b, a / b)
    }

    /// `sqrtsd` (and the `sqrt` call the Release build makes for a negative
    /// argument, which returns the same default NaN).
    pub(crate) fn sqrt(a: f64) -> f64 {
        if a.is_nan() {
            f64::from_bits(a.to_bits() | QUIET)
        } else if a < 0.0 {
            default_nan()
        } else {
            a.sqrt()
        }
    }

    /// `maxsd x, 1.0` as `std::max(1.0, x)` compiles: `x` when `x > 1`,
    /// otherwise `1.0` (so a NaN becomes `1.0`).
    pub(crate) fn max_one(x: f64) -> f64 {
        if x > 1.0 { x } else { 1.0 }
    }

    /// `cvtss2sd`: a float widened to double, a NaN keeping its sign and
    /// payload and becoming quiet.
    pub(crate) fn widen(x: f32) -> f64 {
        if x.is_nan() {
            let bits = x.to_bits();
            let sign = u64::from(bits >> 31) << 63;
            let payload = u64::from(bits & 0x003f_ffff) << 29;
            f64::from_bits(sign | 0x7ff0_0000_0000_0000 | QUIET | payload)
        } else {
            f64::from(x)
        }
    }

    /// `addss`: `a + b` with `a` the first operand.
    pub(crate) fn add32(a: f32, b: f32) -> f32 {
        let result = a + b;
        if a.is_nan() {
            f32::from_bits(a.to_bits() | QUIET_F32)
        } else if b.is_nan() {
            f32::from_bits(b.to_bits() | QUIET_F32)
        } else if result.is_nan() {
            f32::from_bits(DEFAULT_NAN_F32)
        } else {
            result
        }
    }

    /// `divss`: `a / b`.
    pub(crate) fn div32(a: f32, b: f32) -> f32 {
        let result = a / b;
        if a.is_nan() {
            f32::from_bits(a.to_bits() | QUIET_F32)
        } else if b.is_nan() {
            f32::from_bits(b.to_bits() | QUIET_F32)
        } else if result.is_nan() {
            f32::from_bits(DEFAULT_NAN_F32)
        } else {
            result
        }
    }

    /// 32-bit `cvttsd2si`: truncation toward zero, and `i32::MIN` (the integer
    /// indefinite) for NaN and for every value whose truncation does not fit.
    pub(crate) fn cvttsd2si32(x: f64) -> i32 {
        if x > -2_147_483_649.0 && x < 2_147_483_648.0 {
            x as i32
        } else {
            i32::MIN
        }
    }

    /// 64-bit `cvttsd2si`: truncation toward zero, and `i64::MIN` for NaN and
    /// out-of-range values.
    pub(crate) fn cvttsd2si64(x: f64) -> i64 {
        if (-9_223_372_036_854_775_808.0..9_223_372_036_854_775_808.0).contains(&x) {
            x as i64
        } else {
            i64::MIN
        }
    }

    /// GCC's double to `unsigned long` conversion: `comisd 2^63`; below (or
    /// NaN) a signed `cvttsd2si`, otherwise `subsd 2^63; cvttsd2si; btc 63`.
    pub(crate) fn f64_to_u64(x: f64) -> u64 {
        const TWO_63: f64 = 9_223_372_036_854_775_808.0;
        if x >= TWO_63 {
            (cvttsd2si64(x - TWO_63) as u64) ^ (1 << 63)
        } else {
            cvttsd2si64(x) as u64
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn nan_rules_follow_the_first_operand() {
            let positive = f64::from_bits(0x7ff8_0000_0000_0000);
            let negative = f64::from_bits(DEFAULT_NAN);
            assert_eq!(add(negative, positive).to_bits(), DEFAULT_NAN);
            assert_eq!(add(positive, negative).to_bits(), positive.to_bits());
            assert_eq!(add(1.0, positive).to_bits(), positive.to_bits());
            assert_eq!(sub(f64::INFINITY, f64::INFINITY).to_bits(), DEFAULT_NAN);
            assert_eq!(div(0.0, 0.0).to_bits(), DEFAULT_NAN);
            assert_eq!(sqrt(-1.0).to_bits(), DEFAULT_NAN);
            assert_eq!(sqrt(-0.0).to_bits(), (-0.0_f64).to_bits());
            assert_eq!(max_one(f64::NAN), 1.0);
            assert_eq!(div32(0.0, 0.0).to_bits(), DEFAULT_NAN_F32);
            assert_eq!(
                add32(f32::INFINITY, f32::NEG_INFINITY).to_bits(),
                DEFAULT_NAN_F32
            );
            assert_eq!(
                widen(f32::from_bits(0x7fc0_0000)).to_bits(),
                0x7ff8_0000_0000_0000
            );
            assert_eq!(widen(f32::from_bits(0xffc0_0000)).to_bits(), DEFAULT_NAN);
            assert_eq!(
                widen(f32::from_bits(0x7f80_0001)).to_bits(),
                0x7ff8_0000_2000_0000
            );
        }

        #[test]
        fn conversions_return_the_integer_indefinite() {
            assert_eq!(cvttsd2si32(2_147_483_647.9), i32::MAX);
            assert_eq!(cvttsd2si32(2_147_483_648.0), i32::MIN);
            assert_eq!(cvttsd2si32(-2_147_483_648.9), i32::MIN);
            assert_eq!(cvttsd2si32(-2_147_483_649.0), i32::MIN);
            assert_eq!(cvttsd2si32(-0.9), 0);
            assert_eq!(cvttsd2si32(f64::NAN), i32::MIN);
            assert_eq!(cvttsd2si64(f64::NEG_INFINITY), i64::MIN);
            assert_eq!(f64_to_u64(-0.5), 0);
            assert_eq!(f64_to_u64(-3.0), u64::MAX - 2);
            assert_eq!(f64_to_u64(f64::NAN), 1 << 63);
            assert_eq!(f64_to_u64(f64::INFINITY), 0);
            assert_eq!(f64_to_u64(1e30), 0);
            assert_eq!(f64_to_u64(9_223_372_036_854_775_808.0), 1 << 63);
            assert_eq!(f64_to_u64(18_446_744_073_709_549_568.0), u64::MAX - 2047);
        }
    }
}
