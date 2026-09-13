// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A complex FFT and the packed real transform built on it.
//!
//! Native replacement for the two evergreen entry points that
//! `MATH/STATISTICS/KernelDensityEstimation.cpp` reaches for through
//! `Evergreen/evergreen.hpp` and `FFT/FFT.hpp`:
//! `evergreen::real_fft<evergreen::DIF, false, false, true>` and
//! `evergreen::real_ifft<evergreen::DIF, false, false>`. Only those two, plus
//! `evergreen::Tensor`, `evergreen::Vector` and `evergreen::cpx`, are actually
//! called; the belief-propagation machinery the umbrella header pulls in is not
//! used by any ported code. See `docs/FFT_SUPPORT.md`.
//!
//! # The complex transform comes from `rustfft`
//!
//! The C++ takes its FFT from a third-party library, so this port does too
//! rather than reimplementing one: [`crate::math::fft::fft_in_place`] and
//! [`crate::math::fft::ifft_in_place`] delegate to `rustfft` 6.4.1, and
//! [`crate::math::fft::Complex`] is its `num_complex::Complex<f64>`. An earlier
//! version of this module hand-rolled a radix-2 decimation-in-frequency cascade
//! and its own complex type; both are gone.
//!
//! **Always through `FftPlannerScalar`, never the default planner.** `FftPlanner`
//! chooses AVX, SSE or NEON code at runtime. Measured on the same inputs, its
//! output differed from the scalar path in 92-99% of components, with relative
//! differences up to 4.5e-6 at a million points, while buying only 1.0-1.1x on
//! the power-of-two lengths kernel density estimation uses (1.8x at best, on
//! prime lengths). Scalar keeps a kernel density estimate identical on every
//! machine, which is worth far more here than that speed.
//!
//! Everything here is `f64`, as `evergreen::cpx` is a pair of `double`.
//!
//! # What the source computes
//!
//! The forward transform carries no normalisation and the inverse divides by the
//! transform length, so the pair is
//!
//! ```text
//! X[k] = sum_n x[n] exp(-2 pi i k n / N)
//! x[n] = (1/N) sum_k X[k] exp(+2 pi i k n / N)
//! ```
//!
//! which is exactly `rustfft`'s forward transform and its unnormalised inverse
//! scaled by `1/N`.
//!
//! [`crate::math::fft::real_fft`] uses the packing trick `DIF::real_fft1d_packed`
//! uses: a real signal of length `N` is read as `N/2` complex values
//! `x[2j] + i x[2j+1]`, transformed at half length, and unpacked by the
//! equivalent of `RealFFTPostprocessor`. That packing is OpenMS-side behaviour,
//! not library behaviour, so it stays ported here. The result is the first
//! `N/2 + 1` bins; the rest follow from `X[N-k] = conj(X[k])`.
//!
//! # Differences from the source
//!
//! * **Not bit-identical to evergreen, by construction.** evergreen advances a
//!   running twiddle by a recurrence whose error accumulates across a stage;
//!   `rustfft` uses its own planned algorithms. Neither matches evergreen's bit
//!   pattern, and there is no oracle for it, so matching it was not attempted.
//!   [`crate::math::fft::fft`] is pinned against a naive `O(n^2)` DFT instead,
//!   which checks the answer rather than the arithmetic order.
//! * **The complex transform accepts any length.** evergreen derives its length
//!   from `integer_log2(len) = round(log2(len))`, whose `SHAPE_CHECK` assertion
//!   is compiled out of a release build, and dispatches through
//!   `LinearTemplateSearch<0, FFT1D_MAX_LOG_N=16, ...>`, whose terminal case runs
//!   the length-65536 transform regardless. A non-power-of-two length, or one
//!   above 65536, therefore transforms a different number of points than asked
//!   and reports nothing. The earlier hand-rolled kernel refused both. `rustfft`
//!   transforms arbitrary lengths correctly, so [`crate::math::fft::fft`] and
//!   [`crate::math::fft::ifft`] now accept any length from 1 to
//!   [`crate::math::fft::MAX_LEN`]. The packed real transform still requires a
//!   power of two, because its `N/4` unpacking loop does.
//! * **Serial.** Neither `FFT.hpp` nor its callers carry a `#pragma omp`.

use crate::{Error, Result};
use rustfft::FftPlannerScalar;
use std::f64::consts::PI;

/// Largest transform length, `2^24` points.
///
/// Native guard. A transform allocates one buffer of this length; the source
/// allocates whatever it is handed and mis-dispatches above `2^16` (see the
/// module note). Chosen well above the source's reach so that nothing this port
/// accepts is something the source would have computed correctly.
pub const MAX_LEN: usize = 1 << 24;

/// A complex number in `f64`, as `evergreen::cpx` is a pair of `double`.
///
/// This is `rustfft`'s re-export of `num_complex::Complex<f64>`, so values pass
/// to and from the transform without conversion. Construct a real value with
/// `Complex::new(re, 0.0)`; the modulus is [`Complex::norm`].
pub type Complex = rustfft::num_complex::Complex<f64>;

/// `exp(-2 pi i k / n)`, the forward twiddle factor.
///
/// The four quarter-turns are returned exactly rather than through `cos`/`sin`,
/// which would give `6.1e-17` where the answer is zero. The source's running
/// recurrence reaches those points with its own accumulated error.
fn twiddle(k: usize, n: usize) -> Complex {
    // Only the packed real transform calls this, and it requires a power of
    // two, so the quarter-turn test is exact.
    let quarters = 4 * k;
    if quarters % n == 0 {
        return match (quarters / n) % 4 {
            0 => Complex::new(1.0, 0.0),
            1 => Complex::new(0.0, -1.0),
            2 => Complex::new(-1.0, 0.0),
            _ => Complex::new(0.0, 1.0),
        };
    }
    let angle = -2.0 * PI * (k as f64) / (n as f64);
    Complex::new(angle.cos(), angle.sin())
}

/// Reject a length that is zero or above [`MAX_LEN`]. Any other length is a
/// valid complex transform.
fn check_complex_length(len: usize, what: &str) -> Result<()> {
    if len == 0 {
        return Err(Error::InvalidValue(format!("{what} must be positive")));
    }
    if len > MAX_LEN {
        return Err(Error::InvalidRange(format!(
            "{what} {len} exceeds the maximum transform length {MAX_LEN}"
        )));
    }
    Ok(())
}

/// Reject a length the packed real transform cannot unpack: zero, not a power
/// of two, or above [`MAX_LEN`]. Its `N/4` unpacking loop needs a power of two.
fn check_length(len: usize, what: &str) -> Result<()> {
    if len == 0 || !len.is_power_of_two() {
        return Err(Error::InvalidValue(format!(
            "{what} must be a power of two, got {len}"
        )));
    }
    check_complex_length(len, what)
}

/// Forward complex FFT in place, `evergreen::apply_fft<DIF, true, false>`.
///
/// On return `data[k]` holds `sum_n data[n] exp(-2 pi i k n / N)` in natural
/// frequency order. No normalisation is applied, matching the source.
///
/// Computed by `rustfft` through `FftPlannerScalar`, so the result is the same
/// on every CPU; see the module documentation for why the SIMD planner is not
/// used.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the length is zero and
/// [`Error::InvalidRange`] when it exceeds [`MAX_LEN`]. Any other length,
/// including a non-power-of-two, is transformed correctly — unlike the source,
/// which rounds `log2(len)` and silently transforms a different number of points.
pub fn fft_in_place(data: &mut [Complex]) -> Result<()> {
    check_complex_length(data.len(), "FFT length")?;
    // ponytail: a planner per call; cache one per thread if transforms become hot.
    FftPlannerScalar::new()
        .plan_fft_forward(data.len())
        .process(data);
    Ok(())
}

/// Forward complex FFT, returning a new vector.
///
/// # Errors
///
/// As [`fft_in_place`].
pub fn fft(data: &[Complex]) -> Result<Vec<Complex>> {
    let mut out = data.to_vec();
    fft_in_place(&mut out)?;
    Ok(out)
}

/// Inverse complex FFT in place, `evergreen::apply_ifft<DIF, true, false>`.
///
/// On return `data[n]` holds `(1/N) sum_k data[k] exp(+2 pi i k n / N)`, so
/// [`ifft_in_place`] undoes [`fft_in_place`]. `rustfft`'s inverse is
/// unnormalised; the `1/N` scale is applied here, as `NDFFTEnvironment::SingleIFFT1D`
/// applies it in the source.
///
/// # Errors
///
/// As [`fft_in_place`].
pub fn ifft_in_place(data: &mut [Complex]) -> Result<()> {
    check_complex_length(data.len(), "inverse FFT length")?;
    FftPlannerScalar::new()
        .plan_fft_inverse(data.len())
        .process(data);
    let scale = 1.0 / data.len() as f64;
    for value in data.iter_mut() {
        *value *= scale;
    }
    Ok(())
}

/// Inverse complex FFT, returning a new vector.
///
/// # Errors
///
/// As [`fft_in_place`].
pub fn ifft(data: &[Complex]) -> Result<Vec<Complex>> {
    let mut out = data.to_vec();
    ifft_in_place(&mut out)?;
    Ok(out)
}

/// Forward FFT of a real signal, returning its `N/2 + 1` distinct bins.
///
/// Port of `evergreen::real_fft<DIF, false, false, true>` specialised to one
/// dimension, which is the only way `KernelDensityEstimation.cpp` calls it.
/// `values` is read as `N/2` complex numbers `values[2j] + i values[2j+1]`,
/// transformed at half length, and unpacked by the equivalent of
/// `RealFFTPostprocessor::apply`. Bin `k` of the result is
/// `sum_n values[n] exp(-2 pi i k n / N)` for `k` in `0..=N/2`; the bins above
/// `N/2` are `conj` of those below and are not returned. Bins `0` and `N/2`
/// are real, and their imaginary parts are returned as exact zeros.
///
/// The result is **unscaled**, as the source's is.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `values` is empty, has a length that is
/// not a power of two, or contains a non-finite value, and
/// [`Error::InvalidRange`] when it is longer than [`MAX_LEN`]. The source
/// checks none of these; a non-finite input propagates NaN through every bin
/// silently.
pub fn real_fft(values: &[f64]) -> Result<Vec<Complex>> {
    check_length(values.len(), "real FFT length")?;
    if values.iter().any(|v| !v.is_finite()) {
        return Err(Error::InvalidValue(
            "real FFT input must be finite".to_string(),
        ));
    }
    let n = values.len();
    let half = n / 2;
    if n == 1 {
        // A one-point real transform is its own value; the source's
        // `DIF<0>::real_fft1d_packed` is a no-op on the single packed element.
        return Ok(vec![Complex::new(values[0], 0.0)]);
    }
    let mut packed: Vec<Complex> = (0..half)
        .map(|j| Complex::new(values[2 * j], values[2 * j + 1]))
        .collect();
    fft_in_place(&mut packed)?;

    let mut out = vec![Complex::ZERO; half + 1];
    // `RealFFTPostprocessor::apply`, whose first action is to overwrite the
    // zero bin with the two real end bins.
    let bias = packed[0];
    out[0] = Complex::new(bias.re + bias.im, 0.0);
    out[half] = Complex::new(bias.re - bias.im, 0.0);
    for k in 1..=n / 4 {
        let back = packed[half - k].conj();
        let x1 = (packed[k] + back) * 0.5;
        let x2 = (packed[k] - back) * 0.5;
        let w = twiddle(k, n);
        // The source multiplies by `cpx{w.i, -w.r}`, which is `-i w`.
        let temp = x2 * Complex::new(w.im, -w.re);
        out[k] = x1 + temp;
        out[half - k] = (x1 - temp).conj();
    }
    Ok(out)
}

/// Inverse of [`real_fft`]: rebuild a real signal of `length` points from its
/// `length / 2 + 1` distinct bins.
///
/// Port of `evergreen::real_ifft<DIF, false, false>` in one dimension. The
/// imaginary parts of bins `0` and `length / 2` are ignored, as
/// `RealFFTPostprocessor::apply_inverse` ignores them: a real signal has none
/// there, and the source's own `revRt` never supplies them.
///
/// The result is scaled by `1/length` relative to [`real_fft`], so
/// `real_ifft(real_fft(x), x.len())` reproduces `x`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `length` is zero or not a power of two,
/// when `spectrum` does not have exactly `length / 2 + 1` entries, or when it
/// contains a non-finite value, and [`Error::InvalidRange`] when `length`
/// exceeds [`MAX_LEN`].
pub fn real_ifft(spectrum: &[Complex], length: usize) -> Result<Vec<f64>> {
    check_length(length, "inverse real FFT length")?;
    let half = length / 2;
    let expected = half + 1;
    if spectrum.len() != expected {
        return Err(Error::InvalidValue(format!(
            "inverse real FFT needs {expected} bins for {length} points, got {}",
            spectrum.len()
        )));
    }
    if spectrum.iter().any(|v| !v.is_finite()) {
        return Err(Error::InvalidValue(
            "inverse real FFT input must be finite".to_string(),
        ));
    }
    if length == 1 {
        return Ok(vec![spectrum[0].re]);
    }

    // `RealFFTPostprocessor::apply_inverse`, rebuilding the half-length complex
    // spectrum in place over the first `half` entries.
    let mut packed = vec![Complex::ZERO; half];
    let bias = spectrum[0];
    let last = spectrum[half];
    packed[0] = Complex::new((bias.re + last.re) / 2.0, (bias.re - last.re) / 2.0);
    for k in 1..=length / 4 {
        let from_back = spectrum[half - k].conj();
        let x1 = (spectrum[k] + from_back) * 0.5;
        let temp = (spectrum[k] - from_back) * 0.5;
        let w = twiddle(k, length);
        // The source multiplies by `cpx{w.i, w.r}`, the inverse of the forward
        // `-i w`. It stores index `half - k` before index `k` so that the `k`
        // version wins when the two coincide at `k == length / 4`.
        let x2 = temp * Complex::new(w.im, w.re);
        packed[half - k] = (x1 - x2).conj();
        packed[k] = x1 + x2;
    }

    ifft_in_place(&mut packed)?;
    let mut out = vec![0.0; length];
    for (j, value) in packed.iter().enumerate() {
        out[2 * j] = value.re;
        out[2 * j + 1] = value.im;
    }
    Ok(out)
}
