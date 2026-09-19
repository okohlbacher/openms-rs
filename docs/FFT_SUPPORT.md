# FFT: native replacement for the vendored evergreen transform

[`src/math/fft.rs`](../src/math/fft.rs) provides a one-dimensional FFT, backed by
the `rustfft` crate, in place of the vendored `evergreen` library that core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` pulls into
`MATH/STATISTICS/KernelDensityEstimation.cpp` and
`MATH/STATISTICS/MultipleTesting.cpp`.

Tests: [`tests/fft.rs`](../tests/fft.rs). Manifest:
[`tests/data/math_kde_provenance.json`](../tests/data/math_kde_provenance.json).

This module owns no OpenMS header. It exists because the two `.cpp` files above
include `Evergreen/evergreen.hpp` and then, inside `namespace evergreen`,
`FFT/FFT.hpp` — an umbrella that drags in evergreen's whole belief-propagation
machinery. Reading the call sites shows what is actually used:

| evergreen symbol | used for |
| --- | --- |
| `evergreen::cpx` | a pair of `double` |
| `evergreen::Tensor<T>` | a flat 1-D buffer |
| `evergreen::Vector<unsigned long>` | the one-element shape of that buffer |
| `evergreen::DIF` | the decimation-in-frequency 1-D FFT template |
| `evergreen::real_fft<DIF, false, false, true>` | forward real transform, called once in `forRt` |
| `evergreen::real_ifft<DIF, false, false>` | inverse real transform, called once in `revRt` |

Nothing else. The port is therefore a complex FFT and the packed real transform
built on it, not a tensor library.

## Why the transform comes from a crate

The project's policy is that where the C++ takes something from a third-party
library, the port looks for a suitable crate before writing its own. evergreen is
exactly that case. The complex transform is `rustfft` 6.4.1 — pure Rust, MSRV
1.61, and widely used — and `Complex` is its `num_complex::Complex<f64>`. An
earlier version of this module hand-rolled a radix-2 decimation-in-frequency
cascade and a complex type; both were removed.

It is always used through **`FftPlannerScalar`**, never the default `FftPlanner`.
The default picks AVX, SSE or NEON code at runtime. Measured on identical inputs
on the Linux build node:

| length | components differing in the bits | max relative difference | SIMD speedup |
| --- | --- | --- | --- |
| 1,024 | 92% | 1.4e-11 | 1.0x |
| 65,536 | 98% | 9.0e-9 | 1.1x |
| 1,048,576 | 99% | 4.5e-6 | 1.1x |
| 997 (prime) | 99.8% | 1.3e-9 | 1.3x |
| 100,003 (prime) | 99.9% | 2.2e-7 | 1.8x |

The scalar planner avoids runtime selection of CPU-specific SIMD kernels. The
measurements above compare scalar and SIMD execution on the same Linux node;
they do not establish bit-identical results across machines. Cross-machine
bitwise reproducibility remains the intended contract, pending comparisons
across the supported architectures and toolchains. On the measured power-of-two
lengths used by kernel density estimation, SIMD improved speed by at most 1.1x.

## API mapping

| Source construct | Native representation |
| --- | --- |
| `cpx` (`src/openms/extern/evergreen/src/FFT/cpx.hpp`) | `Complex`, a type alias for `rustfft::num_complex::Complex<f64>`; construct a real value with `Complex::new(re, 0.0)`, and the modulus is `norm` |
| `cpx::operator*` | `num_complex`'s `Mul` |
| `cpx::conj` | `num_complex`'s `conj` |
| `DIFButterfly<N>::apply` + `RecursiveShuffle<cpx, LOG_N>::apply` | `fft_in_place` / `fft`, delegated to a transform selected by `rustfft::FftPlannerScalar` |
| `NDFFTEnvironment::SingleIFFT1D::apply` | `ifft_in_place` / `ifft` — conjugate, forward transform, conjugate, scale by `1/N`, in that order |
| `real_fft<DIF, false, false, true>(Tensor<double>)` | `real_fft(&[f64]) -> Result<Vec<Complex>>`, returning the `N/2 + 1` distinct bins |
| `real_ifft<DIF, false, false>(Tensor<cpx>)` | `real_ifft(&[Complex], length) -> Result<Vec<f64>>` |
| `DIF::real_fft1d_packed` + `RealFFTPostprocessor<LOG_N>::apply` | the unpacking half of `real_fft` |
| `RealFFTPostprocessor<LOG_N>::apply_inverse` + `DIF::real_ifft1d_packed` | the repacking half of `real_ifft` |
| `Twiddles<N>::advance`, `Twiddles<N>::delta` | `twiddle(k, n)` inside the packed real transform only, evaluated rather than recurred; the complex transform's twiddles are `rustfft`'s own |
| `integer_log2`, `real_length_to_packed_length`, `packed_length_to_real_length` | folded into `check_length` and the `N/2 + 1` arithmetic |
| `shape_to_log_shape`, `Tensor`, `Vector`, `MatrixTranspose`, `LinearTemplateSearch`, the N-dimensional paths, `apply_fft`, `execute_fft`, `fft_convolve` | **not ported**: no ported call site is multidimensional, and a `Vec<Complex>` is the 1-D tensor |
| `DIT` (decimation in time) | **not ported**: `KernelDensityEstimation.cpp` instantiates `DIF` |
| `FFT1D_MAX_LOG_N = 16` | `MAX_LEN = 2^24`, a real ceiling rather than a dispatch limit |

`Complex` is no longer the crate's own type: `rustfft` already depends on
`num_complex`, and re-exporting its type means values cross into the transform
without conversion.

## Preserved source conventions

- **Unnormalised forward transform, `1/N` on the inverse.** `DIFButterfly`
  applies no scale, and `SingleIFFT1D` divides once by the transform length. The
  pair composes to the identity.
- **The real transform's packing.** A length-`N` real signal is read as `N/2`
  complex values `x[2j] + i x[2j+1]`, transformed at half length, and split by
  the `RealFFTPostprocessor` identities
  `X_k = E_k - i W_N^k O_k`, `X_{N/2-k} = conj(E_k + i W_N^k O_k)`,
  `X_0 = Re Z_0 + Im Z_0`, `X_{N/2} = Re Z_0 - Im Z_0`.
- **The `k == N/4` write order.** The source's forward postprocessor stores
  `data[k]` then `data[N/2-k]` and its inverse stores `data[N/2-k]` then
  `data[k]`, with a comment that the order matters where the two coincide. Both
  orders are reproduced; at that index the two expressions are algebraically the
  same value, and the port asserts as much through the round-trip test.
- **`DIF<0>` is a no-op**, so a one-point transform returns its input rather than
  failing.

## Native differences

- **The complex transform is `rustfft`'s, not evergreen's.** evergreen advances a running
  twiddle by `w += w * delta`, with `delta = (cos(t) - 1, -sin(t))` written as
  `-2 sin^2(t/2)` so the increment stays near zero for large `N`. That reduces a
  recurrence's error without removing it. `rustfft` plans its own algorithms, so
  the results are **not bit-identical** to evergreen's in the last places — nor
  were the earlier hand-rolled kernel's. No oracle for
  evergreen's bit pattern exists in this repository, so matching it was not
  attempted; the port is pinned against a naive `O(n^2)` DFT instead, which
  checks the answer rather than a previous implementation's arithmetic order.
- **The complex transform accepts any length.** evergreen silently transforms
  the wrong number of points for a non-power-of-two length (see the defects
  below), so the earlier hand-rolled kernel refused one. `rustfft` transforms
  arbitrary lengths correctly — prime lengths through Bluestein's algorithm — so
  `fft` and `ifft` now accept any length from 1 to `MAX_LEN`, and refuse only zero
  (`Error::InvalidValue`) or an oversized length (`Error::InvalidRange`). The
  packed real transform still requires a power of two, because its `N/4`
  unpacking loop does.
- **Non-finite inputs are refused** by `real_fft` and `real_ifft`, where the
  source propagates `NaN` through every bin.
- **`real_ifft` takes the output length explicitly** and verifies the bin count
  is exactly `length / 2 + 1`, where the source infers the length from the
  buffer and trusts it.
- **Bins `0` and `N/2` come back with an exactly zero imaginary part**, matching
  `RealFFTPostprocessor::apply`, which writes literal zeros there.
- **Serial.** Neither the vendored FFT nor its two ported call sites carries a
  `#pragma omp`.

## Source defects found

1. **`LinearTemplateSearch` falls through silently above `log2(N) = 16`.**
   `FFT.hpp:7` sets `FFT1D_MAX_LOG_N = 16` and `TemplateSearch.hpp:30`
   terminates the search with `assert(v == MAXIMUM); WORKER<MAXIMUM>::apply(...)`.
   With `NDEBUG` — every release build — the assertion vanishes and a transform
   of any larger length runs the **length-65536** transform over the caller's
   buffer. `gridKdeFFT` derives `M = bit_ceil(max(gridsize, n, 512))`, so a
   sample of more than 65536 points, or a `gridsize` above it, silently produces
   a wrong density with no diagnostic. Reachable from `MultipleTesting::lfdr`,
   whose `gridsize` is a caller parameter.
2. **`integer_log2` rounds instead of rejecting.**
   `integer_log2` (`shape_to_log_shape.hpp:4-14`) computes `round(log2(val))`
   (`:8`) and guards the power-of-two assertion behind `#ifdef SHAPE_CHECK`
   (`:9-11`), which the OpenMS build does not define. A non-power-of-two length therefore transforms a different
   number of points than the caller asked for. `KernelDensityEstimation::forRt`
   documents "rounds up to next power of 2" and implements no such thing, so the
   header and the code disagree as well.

Both are reported in the work package's C++ issue list.

## Checked boundaries and evidence

| Boundary | Where |
| --- | --- |
| `MAX_LEN = 2^24` points per transform | `check_complex_length`, before any allocation |
| complex transform length must be positive | `check_complex_length` |
| packed real transform length must be a power of two | `check_length` |
| `real_ifft` bin count must be `length / 2 + 1` | `real_ifft` |
| every input value must be finite | `real_fft`, `real_ifft` |

Evidence is **tier 4** throughout — there is no evergreen class test in the
OpenMS suite and nothing to transcribe. `tests/fft.rs` asserts:

- the forward transform against a naive `O(n^2)` DFT at every power-of-two
  length from `1` to `1024`, on a deterministic pseudo-random stream;
- the forward transform and its round trip against the naive DFT at lengths the
  earlier kernel refused — 3, 5, 6, 7, 97, 100, 257 and 1000 — including primes,
  which exercise `rustfft`'s Bluestein path;
- a forward/inverse round trip over the same lengths;
- the packed real transform against the naive DFT of the same signal, plus its
  own round trip;
- exact identities: `fft([1,2,3,4]) == [10, -2+2i, -2, -2-2i]`, an impulse
  transforming to a constant, a constant transforming to an impulse of height
  `N` with **exactly** zero elsewhere, `X_0` as the sum and `X_{N/2}` as the
  alternating sum;
- the circular convolution theorem, which is the identity the kernel density
  estimator depends on;
- `fft` and `fft_in_place` agreeing bit for bit;
- every refusal above.

`tests/kernel_density.rs` is the fidelity check that decided the swap: its
upstream class-test expectations pass unchanged on `rustfft`'s scalar planner,
exactly as they did on the hand-rolled kernel.

`tests/kernel_density.rs` adds an indirect check: `for_rt`'s Munro packing is
compared entry by entry against `real_fft`'s output.
