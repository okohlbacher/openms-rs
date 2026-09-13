// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $

//! FFT-based Gaussian kernel density estimation, after Silverman (1982).
//!
//! Port of `src/openms/include/OpenMS/MATH/STATISTICS/KernelDensityEstimation.h`
//! and its translation unit. See `docs/KERNEL_DENSITY_SUPPORT.md`.
//!
//! The estimator bins the sample onto a regular grid, multiplies the grid's
//! transform by the Gaussian kernel's analytic frequency response and transforms
//! back, so the cost is `O(n + M log M)` rather than the `O(n M)` of evaluating
//! a kernel at every grid point. That is the algorithm of
//!
//! > B. W. Silverman (1982). Algorithm AS 176: Kernel density estimation using
//! > the Fast Fourier Transform. J. R. Statist. Soc. C 31(1):93-99.
//!
//! as reimplemented in `statsmodels.nonparametric` and, through it, in
//! PyProphet — whose conventions the source follows where they differ from the
//! paper, and which are called out at the items that carry them.
//!
//! The transforms are [`crate::math::fft`]; the source calls evergreen.
//!
//! # Frequency-domain layout
//!
//! [`crate::math::kernel_density::for_rt`] and
//! [`crate::math::kernel_density::rev_rt`] exchange **Munro-packed** vectors:
//! `M` reals holding `[Re Y_0 .. Re Y_{M/2}, Im Y_1 .. Im Y_{M/2-1}]`. The two
//! bins that are real for a real signal, `0` and `M/2`, contribute no imaginary
//! entry, so `M/2 + 1` real parts and `M/2 - 1` imaginary parts fill exactly
//! `M` slots.
//!
//! # Scaling
//!
//! `statsmodels`' `forrt` divides by `M` and its `revrt` multiplies by `M`; the
//! source does neither, so its `for_rt` is the plain unscaled transform and its
//! `rev_rt` the plain inverse. Both pairs compose to the identity and the
//! kernel product between them is unaffected, which is why the estimates agree.
//! The header's claim that `revRt` is "scaled by multiplying by M" describes
//! `statsmodels`, not the code beneath it, and is recorded as a source defect
//! rather than reproduced — the class test's own round-trip section would fail
//! if it were true.
//!
//! # Differences from the source
//!
//! * Every entry point returns [`crate::Result`] and refuses the degenerate inputs the
//!   source divides by: a zero-width grid, a zero or non-finite `range`, a
//!   non-power-of-two transform length. The per-item docs name each one.
//! * Work is bounded by [`crate::math::kernel_density::MAX_ITEMS`] and
//!   [`crate::math::kernel_density::MAX_GRID`].
//! * Serial, as the source is: `KernelDensityEstimation.cpp` carries no
//!   `#pragma omp`.

use crate::math::fft::{Complex, real_fft, real_ifft};
use crate::math::statistic_functions::quantile;
use crate::{Error, Result};
use std::f64::consts::PI;

/// Maximum number of sample values one call may consume.
///
/// Native guard; the source allocates whatever it is handed.
pub const MAX_ITEMS: usize = 50_000_000;

/// Maximum number of grid points, `2^22`.
///
/// Native guard on the transform length `M` that
/// [`grid_kde_fft`] derives from `gridsize` and the sample
/// size. The source has no ceiling, and above `2^16` its FFT dispatch silently
/// falls back to a length-65536 transform; see [`crate::math::fft`].
pub const MAX_GRID: usize = 1 << 22;

/// The source's default grid size for [`grid_kde_fft`] and
/// [`kde_fft_eval`].
pub const DEFAULT_GRIDSIZE: usize = 512;

/// The source's default grid extension, in bandwidths, beyond the data range.
pub const DEFAULT_CUT: f64 = 3.0;

/// The smallest grid the source will use, irrespective of `gridsize`.
///
/// `gridKdeFFT` raises its target to `max(gridsize, n, 512)` before rounding up
/// to a power of two, so a caller asking for 64 points on 20 samples still gets
/// 512.
pub const MIN_GRID: usize = 512;

/// The divisor in the "nrd0" rule that turns an interquartile range into a
/// comparable standard deviation, `1.34`.
///
/// The exact constant for a normal distribution is `2 Phi^-1(0.75) = 1.349`;
/// R's `bw.nrd0`, `statsmodels` and the source all use the rounded `1.34`, so
/// the port does too.
pub const NRD0_IQR_DIVISOR: f64 = 1.34;

/// `!(a > b)`: true when `a <= b` **and** when either value is `NaN`.
///
/// Spelled through `partial_cmp` because clippy rejects the negated comparison
/// on a partially ordered type. The `NaN` arm is load-bearing at every call
/// site: the source writes `if (!(x > 0.0))` precisely so that a `NaN` takes
/// the guarded branch, and `x <= b` would not.
fn not_greater(a: f64, b: f64) -> bool {
    !matches!(a.partial_cmp(&b), Some(std::cmp::Ordering::Greater))
}

fn bad(message: &str) -> Error {
    Error::InvalidValue(message.to_string())
}

fn check_items(len: usize) -> Result<()> {
    if len > MAX_ITEMS {
        return Err(Error::InvalidRange(format!(
            "kernel density input of {len} values exceeds the maximum {MAX_ITEMS}"
        )));
    }
    Ok(())
}

/// Bandwidth by Silverman's "nrd0" rule of thumb.
///
/// Computes `0.9 * min(sd, IQR / 1.34) * n^(-1/5)` over the finite values of
/// `x`, where `sd` is the sample standard deviation with one degree of freedom
/// removed and the quartiles are the interpolating
/// [`crate::math::statistic_functions::quantile`] at `0.25` and `0.75` — the
/// `numpy.percentile` convention, not the median-of-halves convention of
/// [`crate::math::statistic_functions::quantile1st_sorted`].
///
/// Non-finite values are dropped before anything is computed, and fewer than
/// two survivors yield `0.0` rather than an error: the source returns `0.0`
/// there and its callers treat a zero bandwidth as "no estimate".
///
/// Equivalent to R's `bw.nrd0()` and `statsmodels`' `bw_silverman()`.
///
/// # The fallback chain
///
/// When `min(sd, IQR / 1.34)` is not strictly positive — a constant sample, or
/// one whose middle half is constant — the source falls back in PyProphet's
/// order: the standard deviation, then the absolute value of the **smallest**
/// value (the sample has been sorted by that point, so PyProphet's `x[0]` is
/// the minimum, not the first input), then `1.0`. A sample of ten zeros
/// therefore yields `0.9 * 10^(-1/5)` rather than a zero-width grid downstream.
///
/// Reference: B. W. Silverman (1986). Density Estimation for Statistics and
/// Data Analysis. Chapman & Hall/CRC.
///
/// # Errors
///
/// Returns [`Error::InvalidRange`] when `x` is longer than
/// [`MAX_ITEMS`]. The source has no ceiling.
pub fn bw_nrd0(x: &[f64]) -> Result<f64> {
    check_items(x.len())?;
    let mut finite: Vec<f64> = x.iter().copied().filter(|v| v.is_finite()).collect();
    let n = finite.len();
    if n < 2 {
        return Ok(0.0);
    }
    // Sample standard deviation, accumulated in input order as the source does.
    let mut mean = 0.0;
    for value in &finite {
        mean += *value;
    }
    mean /= n as f64;
    let mut var = 0.0;
    for value in &finite {
        let d = *value - mean;
        var += d * d;
    }
    var /= (n - 1) as f64;
    let sd = var.sqrt();

    finite.sort_by(f64::total_cmp);
    let q25 = quantile(&finite, 0.25)?;
    let q75 = quantile(&finite, 0.75)?;
    let iqr = q75 - q25;

    let scaled_iqr = iqr / NRD0_IQR_DIVISOR;
    let mut lo = if scaled_iqr < sd { scaled_iqr } else { sd };
    if not_greater(lo, 0.0) {
        if sd > 0.0 {
            lo = sd;
        } else {
            lo = finite[0].abs();
        }
    }
    if not_greater(lo, 0.0) {
        lo = 1.0;
    }
    Ok(0.9 * lo * (n as f64).powf(-0.2))
}

/// Count `x` into `nbins` equal-width bins spanning `[xmin, xmax]`.
///
/// **This is a histogram, not linear binning**, despite the name and the
/// header's description of proportional allocation to the two nearest grid
/// points. Each value lands whole in the bin `floor((v - xmin) / width)` with
/// `width = (xmax - xmin) / nbins`, and `xmax` itself is folded into the last
/// bin. The class test pins that behaviour: `linBin([0.1, 0.4, 0.9], 0, 1, 5)`
/// is `[1, 0, 1, 0, 1]`, which proportional allocation could not produce. True
/// linear binning does exist in the source, in the file-static `fast_linbin`
/// that [`grid_kde_fft`] uses, and is reproduced there.
/// The mismatch is recorded in `OpenMS_CPP_ISSUES.md`.
///
/// Values outside `[xmin, xmax]`, and non-finite values, are ignored.
///
/// # Arguments
///
/// * `weights` — one weight per value, or `None` for unit weights. A slice
///   whose length differs from `x` is **ignored** and unit weights are used, as
///   the source's `weights != nullptr && weights->size() == x.size()` test does;
///   the port keeps that rather than erroring, because callers rely on passing a
///   possibly-empty vector.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `nbins` is zero, when `xmax` is not
/// greater than `xmin`, or when either bound is not finite — the first two are
/// the source's `std::invalid_argument`, the third is native — and
/// [`Error::InvalidRange`] when `x` is longer than [`MAX_ITEMS`]
/// or `nbins` exceeds [`MAX_GRID`].
pub fn lin_bin(
    x: &[f64],
    xmin: f64,
    xmax: f64,
    nbins: usize,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>> {
    if nbins == 0 {
        return Err(bad("linBin: nbins must be > 0"));
    }
    if nbins > MAX_GRID {
        return Err(Error::InvalidRange(format!(
            "linBin: {nbins} bins exceeds the maximum {MAX_GRID}"
        )));
    }
    check_items(x.len())?;
    let mut bins = vec![0.0; nbins];
    if x.is_empty() {
        // The source returns the zeroed bins before validating the bounds.
        return Ok(bins);
    }
    if !xmin.is_finite() || !xmax.is_finite() {
        return Err(bad("linBin: bounds must be finite"));
    }
    if not_greater(xmax, xmin) {
        return Err(bad("linBin: xmax must be > xmin"));
    }
    let width = (xmax - xmin) / nbins as f64;
    if not_greater(width, 0.0) {
        // Reachable only when the span underflows; the source returns zeros.
        return Ok(bins);
    }
    let use_weights = weights.is_some_and(|w| w.len() == x.len());
    for (i, value) in x.iter().copied().enumerate() {
        if !value.is_finite() || value < xmin || value > xmax {
            continue;
        }
        let position = ((value - xmin) / width).floor();
        // `position` is in [0, nbins] because `value` is within the bounds.
        let mut index = position as usize;
        if index >= nbins {
            index = nbins - 1;
        }
        let w = if use_weights {
            // `use_weights` proved the lengths equal.
            weights.map_or(1.0, |values| values[i])
        } else {
            1.0
        };
        bins[index] += w;
    }
    Ok(bins)
}

/// True linear binning onto `M` points of `linspace(a, b, M)`, the source's
/// file-static `fast_linbin`.
///
/// Each value is split between its two neighbouring grid points in proportion
/// to its distance from them, which is what makes the FFT estimate smooth in
/// the sample rather than in the bin edges. Values outside `[a, b]` and
/// non-finite values are dropped, and a value at or past the last interval goes
/// whole into the last point.
fn fast_lin_bin(x: &[f64], a: f64, b: f64, m: usize) -> Vec<f64> {
    let mut bins = vec![0.0; m];
    if x.is_empty() || not_greater(b, a) || m < 2 {
        return bins;
    }
    // `linspace(a, b, m)` has m - 1 intervals, not m.
    let delta = (b - a) / (m - 1) as f64;
    for value in x.iter().copied() {
        if !value.is_finite() || value < a || value > b {
            continue;
        }
        let position = (value - a) / delta;
        let left_edge = position.floor();
        if left_edge.is_nan() || left_edge < 0.0 || left_edge >= (m - 1) as f64 {
            bins[m - 1] += 1.0;
            continue;
        }
        let left = left_edge as usize;
        let frac = position - left_edge;
        bins[left] += 1.0 - frac;
        bins[left + 1] += frac;
    }
    bins
}

/// Forward real FFT of `values`, zero-padded or truncated to `m` points, in
/// Munro-packed form.
///
/// The result is unscaled. Non-finite entries of `values` are read as zero,
/// which is what the source's `std::isfinite` test in the padding loop does.
///
/// # Arguments
///
/// * `m` — transform length; `0` means `values.len()`. Must be a power of two.
///   The header says a length that is not one is "rounded up to the next power
///   of 2"; the implementation does no such thing and hands the raw length to
///   evergreen, whose `integer_log2` rounds `log2` to the nearest integer with
///   its shape assertion compiled out. This port rejects instead. Recorded in
///   `OpenMS_CPP_ISSUES.md`.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when the effective length is not a power of
/// two and [`Error::InvalidRange`] when it exceeds
/// [`crate::math::fft::MAX_LEN`]. An effective length of zero yields an empty
/// vector, as the source does.
pub fn for_rt(values: &[f64], m: usize) -> Result<Vec<f64>> {
    let m = if m == 0 { values.len() } else { m };
    if m == 0 {
        return Ok(Vec::new());
    }
    check_items(m)?;
    let mut padded = vec![0.0; m];
    for (slot, value) in padded.iter_mut().zip(values.iter().copied()) {
        if value.is_finite() {
            *slot = value;
        }
    }
    let spectrum = real_fft(&padded)?;

    let half = m / 2;
    let mut out = vec![0.0; m];
    for (k, bin) in spectrum.iter().enumerate() {
        out[k] = bin.re;
    }
    // `Im Y_0` and `Im Y_{M/2}` are zero for a real signal and are not stored.
    for k in 1..spectrum.len().saturating_sub(1) {
        out[half + k] = spectrum[k].im;
    }
    Ok(out)
}

/// Inverse of [`for_rt`]: rebuild `m` real points from a
/// Munro-packed spectrum.
///
/// # Arguments
///
/// * `m` — output length; `0` means `packed.len()`. Must be a power of two.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `packed` is shorter than `m` — the
/// source's `std::invalid_argument` — when the effective length is not a power
/// of two, or when `packed` contains a non-finite value, and
/// [`Error::InvalidRange`] when the length exceeds
/// [`crate::math::fft::MAX_LEN`]. An effective length of zero yields an empty
/// vector.
pub fn rev_rt(packed: &[f64], m: usize) -> Result<Vec<f64>> {
    let m = if m == 0 { packed.len() } else { m };
    if m == 0 {
        return Ok(Vec::new());
    }
    check_items(m)?;
    if packed.len() < m {
        return Err(bad("revRt: input length must equal M"));
    }
    let half = m / 2;
    let bins = half + 1;
    let mut spectrum = vec![Complex::ZERO; bins];
    for (k, slot) in spectrum.iter_mut().enumerate() {
        let im = if k > 0 && k < bins - 1 {
            packed[half + k]
        } else {
            0.0
        };
        *slot = Complex::new(packed[k], im);
    }
    real_ifft(&spectrum, m)
}

/// The Gaussian kernel's analytic frequency response, in Munro-packed form.
///
/// Silverman's transform: bin `j` is `exp(-2 (pi bw j / range)^2)` divided by
/// the correction `1 - (j pi / M)^2 / 3`, which compensates the linear binning.
/// Computing it in closed form rather than transforming a sampled Gaussian
/// saves one FFT and avoids the aliasing a truncated kernel would introduce.
/// Bin `0` is exactly `1.0`, so the kernel preserves total mass.
///
/// The correction can go non-positive for `j` near `M/2`; the source replaces
/// it with the smallest positive normal `double`, which drives that bin's
/// response to the largest finite value instead of to a sign flip. That
/// substitution is preserved.
///
/// The mirrored second half repeats bins `1..M/2 - 1`, because the kernel is
/// real and even and its packed imaginary slots multiply the signal's.
///
/// # Arguments
///
/// * `bw` — kernel standard deviation; must be finite.
/// * `m` — grid length. `0` yields an empty vector, as in the source.
/// * `range` — width of the spatial domain the kernel is applied over; must be
///   finite and non-zero.
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `bw` is not finite or `range` is not
/// finite or is zero, and [`Error::InvalidRange`] when `m` exceeds
/// [`MAX_GRID`]. The source divides by `range` unchecked and
/// would return an all-`NaN` kernel for `bw == range == 0`.
pub fn silverman_kernel_fft(bw: f64, m: usize, range: f64) -> Result<Vec<f64>> {
    if m == 0 {
        return Ok(Vec::new());
    }
    if m > MAX_GRID {
        return Err(Error::InvalidRange(format!(
            "silvermanKernelFFT: {m} grid points exceeds the maximum {MAX_GRID}"
        )));
    }
    if !bw.is_finite() {
        return Err(bad("silvermanKernelFFT: bandwidth must be finite"));
    }
    if !range.is_finite() || range == 0.0 {
        return Err(bad("silvermanKernelFFT: range must be finite and non-zero"));
    }
    let half = m / 2;
    let fac1 = 2.0 * (PI * bw / range).powi(2);
    let mut factors = vec![0.0; half + 1];
    for (k, slot) in factors.iter_mut().enumerate() {
        let j = k as f64;
        let jfac = j * j * fac1;
        let scaled = j * (PI / m as f64);
        let mut bc = 1.0 - scaled * scaled / 3.0;
        if not_greater(bc, 0.0) {
            bc = f64::MIN_POSITIVE;
        }
        *slot = (-jfac).exp() / bc;
    }
    let mut out = vec![0.0; m];
    out[..=half].copy_from_slice(&factors);
    if half >= 2 {
        out[half + 1..half + half].copy_from_slice(&factors[1..half]);
    }
    Ok(out)
}

/// A density estimate sampled on a regular grid, the source's
/// `std::pair<std::vector<double>, std::vector<double>>`.
///
/// Named fields replace the pair because `first` and `second` do not say which
/// is which; the source's `first` is [`GridKde::density`] and its `second` is
/// [`GridKde::grid`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct GridKde {
    /// Density at each grid point, renormalised so the trapezoid-free sum
    /// `sum(density) * spacing` is one.
    pub density: Vec<f64>,
    /// The grid points themselves, ascending and equally spaced.
    pub grid: Vec<f64>,
}

/// Gaussian kernel density estimate on a regular grid, by FFT convolution.
///
/// The grid spans `min(x) - cut * bw` to `max(x) + cut * bw` in `M` equally
/// spaced points, where `M` is the next power of two at or above
/// `max(gridsize, n, 512)`; a caller's `gridsize` is therefore a lower bound,
/// not the answer. An empty sample centres the grid on zero.
///
/// The sample is linearly binned onto the grid, divided by `spacing * n`,
/// transformed, multiplied by [`silverman_kernel_fft`], transformed back and
/// finally rescaled so the estimate integrates to one. That last rescaling is
/// the source's, and it is why the density may be pulled slightly away from the
/// unrenormalised convolution.
///
/// # Arguments
///
/// * `bw` — kernel bandwidth, typically from [`bw_nrd0`].
/// * `gridsize` — lower bound on the grid length; the source's default is
///   [`DEFAULT_GRIDSIZE`].
/// * `cut` — grid extension in bandwidths beyond the data range; the source's
///   default is [`DEFAULT_CUT`].
///
/// # Errors
///
/// Returns [`Error::InvalidValue`] when `x` holds a non-finite value, when `bw`
/// or `cut` is not finite, or when the resulting grid has zero width — all
/// three are cases the source carries into `std::minmax_element`, a division by
/// zero, or an all-`NaN` kernel without comment. Returns
/// [`Error::InvalidRange`] when `x` exceeds [`MAX_ITEMS`] or the
/// derived grid length exceeds [`MAX_GRID`].
///
/// The source additionally throws `std::runtime_error` on an internal size
/// mismatch between the kernel and the transform; that cannot arise here,
/// because both lengths are the same `M`.
pub fn grid_kde_fft(x: &[f64], bw: f64, gridsize: usize, cut: f64) -> Result<GridKde> {
    check_items(x.len())?;
    if x.iter().any(|v| !v.is_finite()) {
        return Err(bad("gridKdeFFT: sample values must be finite"));
    }
    if !bw.is_finite() || !cut.is_finite() {
        return Err(bad("gridKdeFFT: bandwidth and cut must be finite"));
    }
    let n = x.len();
    let target = gridsize.max(n.max(MIN_GRID));
    if target > MAX_GRID {
        return Err(Error::InvalidRange(format!(
            "gridKdeFFT: a grid of at least {target} points exceeds the maximum {MAX_GRID}"
        )));
    }
    // `target <= MAX_GRID`, itself a power of two, so this cannot overflow.
    let m = target.next_power_of_two();

    let (a, b) = if n > 0 {
        let mut lo = x[0];
        let mut hi = x[0];
        for value in &x[1..] {
            if *value < lo {
                lo = *value;
            }
            if *value > hi {
                hi = *value;
            }
        }
        (lo - cut * bw, hi + cut * bw)
    } else {
        (-cut * bw, cut * bw)
    };
    if !a.is_finite() || !b.is_finite() || not_greater(b, a) {
        return Err(bad(
            "gridKdeFFT: the grid has zero width; check the bandwidth and cut",
        ));
    }
    let delta = (b - a) / (m - 1) as f64;
    if not_greater(delta, 0.0) {
        return Err(bad("gridKdeFFT: the grid spacing underflows"));
    }
    let grid: Vec<f64> = (0..m).map(|i| a + i as f64 * delta).collect();
    let range = b - a;

    let mut binned = fast_lin_bin(x, a, b, m);
    if n > 0 {
        let scale = delta * n as f64;
        for value in binned.iter_mut() {
            *value /= scale;
        }
    }

    let forward = for_rt(&binned, m)?;
    let kernel = silverman_kernel_fft(bw, m, range)?;
    let product: Vec<f64> = kernel
        .iter()
        .zip(forward.iter())
        .map(|(k, y)| k * y)
        .collect();
    let mut density = rev_rt(&product, m)?;

    let mut sum = 0.0;
    for value in &density {
        sum += *value;
    }
    if sum > 0.0 && delta > 0.0 {
        let scale = 1.0 / (sum * delta);
        for value in density.iter_mut() {
            *value *= scale;
        }
    }
    Ok(GridKde { density, grid })
}

/// Kernel density estimate evaluated at the sample points themselves.
///
/// Runs [`grid_kde_fft`] and interpolates the grid estimate to each point of
/// `x` with a natural cubic spline, which is `O(n log M)` rather than the
/// `O(n^2)` of summing kernels. Negative values, which the spline can produce
/// where the estimate is near zero, are clamped to zero as the source clamps
/// them.
///
/// The spline is the same construction as the crate's own port of
/// `MATH/MISC/CubicSpline2d.h` in `src/processing/spline/cubic.rs`, which is
/// the class the source uses; it is rebuilt here as a private helper because
/// the module dependency graph this crate ratchets does not admit an edge from
/// `math` to `processing`. `tests/kernel_density.rs` asserts the two agree on
/// the actual KDE grid, so the duplication is checked rather than assumed.
///
/// # Errors
///
/// As [`grid_kde_fft`], and additionally [`Error::InvalidValue`] when a point
/// of `x` falls outside the grid — which the source reports as
/// `Exception::IllegalArgument` from `CubicSpline2d::eval` — or when the grid
/// estimate is not usable as spline ordinates.
pub fn kde_fft_eval(x: &[f64], bw: f64, gridsize: usize, cut: f64) -> Result<Vec<f64>> {
    let estimate = grid_kde_fft(x, bw, gridsize, cut)?;
    if x.is_empty() {
        return Ok(Vec::new());
    }
    let spline = NaturalCubicSpline::new(&estimate.grid, &estimate.density)?;
    let mut out = Vec::with_capacity(x.len());
    for point in x.iter().copied() {
        let value = spline.eval(point)?;
        out.push(if value < 0.0 { 0.0 } else { value });
    }
    Ok(out)
}

/// The natural cubic spline of `MATH/MISC/CubicSpline2d.h`, private to this
/// module.
///
/// Reproduces `CubicSpline2d::init_` term by term, including its Thomas-style
/// forward sweep and the Horner evaluation `((d t + c) t + b) t + a`.
struct NaturalCubicSpline {
    x: Vec<f64>,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>,
}

impl NaturalCubicSpline {
    fn new(x: &[f64], y: &[f64]) -> Result<Self> {
        if x.len() != y.len() || x.len() < 2 {
            return Err(bad("cubic spline needs matching arrays of 2 or more knots"));
        }
        if x.iter().chain(y.iter()).any(|v| !v.is_finite()) {
            return Err(bad("cubic spline knots must be finite"));
        }
        if x.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(bad("cubic spline abscissae must be strictly increasing"));
        }
        let n = x.len() - 1;
        let h: Vec<f64> = x.windows(2).map(|pair| pair[1] - pair[0]).collect();
        let mut mu = vec![0.0; n];
        let mut z = vec![0.0; n];
        for i in 1..n {
            let span = x[i + 1] - x[i - 1];
            let l = 2.0 * span - h[i - 1] * mu[i - 1];
            mu[i] = h[i] / l;
            z[i] = (3.0 * (y[i + 1] * h[i - 1] - y[i] * span + y[i - 1] * h[i])
                / (h[i - 1] * h[i])
                - h[i - 1] * z[i - 1])
                / l;
        }
        let mut b = vec![0.0; n];
        let mut c = vec![0.0; n + 1];
        let mut d = vec![0.0; n];
        for j in (0..n).rev() {
            c[j] = z[j] - mu[j] * c[j + 1];
            b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
            d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
        }
        if b.iter()
            .chain(c.iter())
            .chain(d.iter())
            .any(|v| !v.is_finite())
        {
            return Err(bad("cubic spline coefficients are not finite"));
        }
        Ok(Self {
            x: x.to_vec(),
            a: y[..n].to_vec(),
            b,
            c,
            d,
        })
    }

    fn eval(&self, x: f64) -> Result<f64> {
        if !x.is_finite() {
            return Err(bad("cubic spline argument must be finite"));
        }
        let first = self.x[0];
        let last = self.x[self.x.len() - 1];
        if x < first || x > last {
            return Err(bad("cubic spline argument is out of range"));
        }
        // `lower_bound`, then step back when the hit is above `x` or is the
        // final knot, whose segment starts one place earlier.
        let mut i = self.x.partition_point(|knot| *knot < x);
        if i >= self.x.len() || self.x[i] > x || x == last {
            // `x >= first` and the abscissae are strictly increasing, so this
            // branch is never reached with `i == 0`; the saturating form is
            // there so a future edit cannot turn it into a panic.
            i = i.saturating_sub(1);
        }
        let dx = x - self.x[i];
        Ok(((self.d[i] * dx + self.c[i]) * dx + self.b[i]) * dx + self.a[i])
    }
}
