// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! A radix-2 decimation-in-frequency complex FFT and the packed real transform
//! built on it.
//!
//! Native replacement for the two evergreen entry points that
//! `MATH/STATISTICS/KernelDensityEstimation.cpp` reaches for through
//! `Evergreen/evergreen.hpp` and `FFT/FFT.hpp`:
//! `evergreen::real_fft<evergreen::DIF, false, false, true>` and
//! `evergreen::real_ifft<evergreen::DIF, false, false>`. Only those two, plus
//! `evergreen::Tensor`, `evergreen::Vector` and `evergreen::cpx`, are actually
//! called; the belief-propagation machinery the umbrella header pulls in is not
//! used by any ported code, so this module is a one-dimensional FFT and nothing
//! more. See `docs/FFT_SUPPORT.md`.
//!
//! Everything here is `f64`, as `evergreen::cpx` is a pair of `double`.
//!
//! # What the source computes
//!
//! `DIF<LOG_N, SHUFFLE>::fft1d` applies `DIFButterfly<N>` — combine-then-twiddle
//! butterflies over halves, recursing into each half — and then a bit-reversal
//! shuffle, which is the textbook decimation-in-frequency arrangement. The
//! forward transform carries no normalisation and the inverse divides by the
//! transform length, so the pair is
//!
//! ```text
//! X[k] = sum_n x[n] exp(-2 pi i k n / N)
//! x[n] = (1/N) sum_k X[k] exp(+2 pi i k n / N)
//! ```
//!
//! [`crate::math::fft::real_fft`] uses the standard packing trick that
//! `DIF::real_fft1d_packed` uses: a real signal of length `N` is read as `N/2`
//! complex values `x[2j] + i x[2j+1]`, transformed at half length, and unpacked
//! by `RealFFTPostprocessor`. The result is the first `N/2 + 1` bins; the rest
//! follow from `X[N-k] = conj(X[k])`.
//!
//! # Differences from the source
//!
//! * **Twiddle factors are evaluated, not recurred.** evergreen advances a
//!   running twiddle by `w += w * delta` with `delta = (cos(t) - 1, -sin(t))`,
//!   a recurrence chosen to keep the increment near zero for large `N`. It is
//!   still a recurrence, and its error accumulates across a stage. This port
//!   evaluates `cos`/`sin` per butterfly and returns the exact values at the
//!   four quarter-turns, which is more accurate and not bit-identical in the
//!   last places. There is no oracle for evergreen's bit pattern, so matching
//!   it was not attempted; [`crate::math::fft::fft`] is instead pinned against
//!   a naive `O(n^2)` DFT, which checks the answer rather than the arithmetic
//!   order.
//! * **Lengths are checked.** evergreen derives its transform length from
//!   `integer_log2(len) = round(log2(len))`, whose `SHAPE_CHECK` assertion is
//!   compiled out of a release build, and dispatches through
//!   `LinearTemplateSearch<0, FFT1D_MAX_LOG_N=16, ...>`, whose terminal case
//!   asserts `v == 16` and then runs the length-65536 transform regardless. A
//!   non-power-of-two length, or one above 65536, therefore transforms a
//!   different number of points than the caller asked for and reports nothing.
//!   This port returns [`crate::Error::InvalidValue`] in both cases; the ceiling is
//!   [`crate::math::fft::MAX_LEN`], above the source's so that a caller is not
//!   silently more restricted, and the failure mode is an error either way.
//! * **Serial.** Neither `FFT.hpp` nor its callers carry a `#pragma omp`.

use crate::{Error, Result};
use std::f64::consts::PI;
use std::ops::{Add, Mul, Sub};

/// Largest transform length, `2^24` points.
///
/// Native guard. A transform allocates one buffer of this length; the source
/// allocates whatever it is handed and mis-dispatches above `2^16` (see the
/// module note). Chosen well above the source's reach so that nothing this port
/// accepts is something the source would have computed correctly.
pub const MAX_LEN: usize = 1 << 24;

/// A complex number, the port of `evergreen::cpx`.
///
/// The source names the parts `r` and `i`; they are spelled out here because a
/// one-letter `i` next to an index `i` reads badly.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part, the source's `cpx::r`.
    pub re: f64,
    /// Imaginary part, the source's `cpx::i`.
    pub im: f64,
}

impl Complex {
    /// The additive identity.
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };

    /// A complex number from its real and imaginary parts.
    pub fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    /// A real number as a complex one, with a zero imaginary part.
    pub fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }

    /// The complex conjugate, the source's `cpx::conj`.
    pub fn conj(self) -> Self {
        Self {
            re: self.re,
            im: -self.im,
        }
    }

    /// The modulus `sqrt(re^2 + im^2)`.
    ///
    /// Uses [`f64::hypot`], which does not overflow where the naive form would.
    pub fn modulus(self) -> f64 {
        self.re.hypot(self.im)
    }

    /// Whether both parts are finite.
    pub fn is_finite(self) -> bool {
        self.re.is_finite() && self.im.is_finite()
    }
}

impl Add for Complex {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl Sub for Complex {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl Mul for Complex {
    type Output = Self;
    /// The source's `cpx::operator*`: the four-multiply schoolbook product, not
    /// Karatsuba's three-multiply form, which rounds differently.
    fn mul(self, rhs: Self) -> Self {
        Self {
            re: self.re * rhs.re - self.im * rhs.im,
            im: self.re * rhs.im + self.im * rhs.re,
        }
    }
}

impl Mul<f64> for Complex {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        Self {
            re: self.re * rhs,
            im: self.im * rhs,
        }
    }
}

/// `exp(-2 pi i k / n)`, the forward twiddle factor.
///
/// The four quarter-turns are returned exactly rather than through `cos`/`sin`,
/// which would give `6.1e-17` where the answer is zero. The source's running
/// recurrence reaches those points with its own accumulated error.
fn twiddle(k: usize, n: usize) -> Complex {
    // n is a power of two here, so the quarter-turn test is exact.
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

/// Reject a length that is not a positive power of two within [`MAX_LEN`].
fn check_length(len: usize, what: &str) -> Result<()> {
    if len == 0 || !len.is_power_of_two() {
        return Err(Error::InvalidValue(format!(
            "{what} must be a power of two, got {len}"
        )));
    }
    if len > MAX_LEN {
        return Err(Error::InvalidRange(format!(
            "{what} {len} exceeds the maximum transform length {MAX_LEN}"
        )));
    }
    Ok(())
}

/// The bit-reversal permutation `RecursiveShuffle` applies after the
/// decimation-in-frequency butterflies.
fn bit_reverse(data: &mut [Complex]) {
    let n = data.len();
    if n < 4 {
        return;
    }
    let shift = usize::BITS - n.trailing_zeros();
    for i in 0..n {
        let j = i.reverse_bits() >> shift;
        if j > i {
            data.swap(i, j);
        }
    }
}

/// The decimation-in-frequency butterfly cascade, without the final shuffle.
///
/// This is `DIFButterfly<N>::apply` flattened from evergreen's recursion into
/// the equivalent loop: the outermost stage spans the whole range, each
/// following stage halves it, and the two recursive calls become the blocks of
/// the next stage. The per-butterfly operation order — sum first, then
/// difference times twiddle — is the source's.
fn dif_butterflies(data: &mut [Complex]) {
    let n = data.len();
    let mut len = n;
    while len >= 2 {
        let half = len / 2;
        let mut start = 0;
        while start < n {
            for i in 0..half {
                let a = data[start + i];
                let b = data[start + i + half];
                data[start + i] = a + b;
                data[start + i + half] = (a - b) * twiddle(i, len);
            }
            start += len;
        }
        len = half;
    }
}

/// Forward complex FFT in place, `evergreen::apply_fft<DIF, true, false>`.
///
/// On return `data[k]` holds `sum_n data[n] exp(-2 pi i k n / N)` in natural
/// frequency order. No normalisation is applied, matching the source.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the length is zero or not a power of
/// two and [`Error::InvalidRange`] when it exceeds [`MAX_LEN`]. The source
/// rounds `log2(len)` and dispatches on the rounded value, so it neither
/// rejects nor reports either case.
pub fn fft_in_place(data: &mut [Complex]) -> Result<()> {
    check_length(data.len(), "FFT length")?;
    dif_butterflies(data);
    bit_reverse(data);
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
/// Conjugate, forward-transform, conjugate and scale by `1/N`, which is exactly
/// what `NDFFTEnvironment::SingleIFFT1D` does and in that order.
///
/// # Errors
///
/// As [`fft_in_place`].
pub fn ifft_in_place(data: &mut [Complex]) -> Result<()> {
    check_length(data.len(), "inverse FFT length")?;
    for value in data.iter_mut() {
        *value = value.conj();
    }
    dif_butterflies(data);
    bit_reverse(data);
    let scale = 1.0 / data.len() as f64;
    for value in data.iter_mut() {
        *value = value.conj() * scale;
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
        return Ok(vec![Complex::real(values[0])]);
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
