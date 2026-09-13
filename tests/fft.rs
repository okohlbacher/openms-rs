// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//! Tests for `src/math/fft.rs`, the native replacement for the two evergreen
//! FFT entry points `MATH/STATISTICS/KernelDensityEstimation.cpp` calls.
//!
//! evergreen has no class test in the OpenMS suite — it is a vendored library
//! under `src/openms/extern/evergreen` — so there are no literals to
//! transcribe. Every expectation here is tier 4: derived from the definition of
//! the discrete Fourier transform, checked against a naive `O(n^2)` DFT written
//! out below, or asserted as an exact algebraic identity (`X_0` is the sum,
//! `X_{N/2}` is the alternating sum, a delta transforms to a constant, a
//! constant transforms to a delta). That is a stronger check than a transcribed
//! value would be: it tests the answer rather than a previous implementation's
//! arithmetic order.

use openms::Error;
use openms::math::fft::{Complex, MAX_LEN, fft, fft_in_place, ifft, real_fft, real_ifft};
use std::f64::consts::PI;

/// The definition, `X_k = sum_n x_n exp(-2 pi i k n / N)`, evaluated directly.
fn naive_dft(data: &[Complex]) -> Vec<Complex> {
    let n = data.len();
    (0..n)
        .map(|k| {
            let mut acc = Complex::ZERO;
            for (index, value) in data.iter().enumerate() {
                let angle = -2.0 * PI * (k as f64) * (index as f64) / (n as f64);
                let (sin, cos) = angle.sin_cos();
                acc += *value * Complex::new(cos, sin);
            }
            acc
        })
        .collect()
}

/// A deterministic, reproducible pseudo-random stream. `xorshift64*`, chosen so
/// the test needs no dependency and always runs on the same numbers.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// A value in `[-1, 1)`.
    fn next_f64(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let bits = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D);
        ((bits >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

fn assert_close(actual: Complex, expected: Complex, scale: f64, what: &str) {
    let tolerance = 1e-12 * scale.max(1.0);
    assert!(
        (actual.re - expected.re).abs() <= tolerance
            && (actual.im - expected.im).abs() <= tolerance,
        "{what}: ({}, {}) != ({}, {})",
        actual.re,
        actual.im,
        expected.re,
        expected.im
    );
}

#[test]
fn forward_transform_matches_the_naive_dft() {
    let mut rng = Rng::new(0x9E37_79B9_7F4A_7C15);
    for log_n in 0..=10u32 {
        let n = 1usize << log_n;
        let input: Vec<Complex> = (0..n)
            .map(|_| Complex::new(rng.next_f64(), rng.next_f64()))
            .collect();
        let fast = fft(&input).expect("power-of-two length");
        let slow = naive_dft(&input);
        for (k, (got, want)) in fast.iter().zip(slow.iter()).enumerate() {
            assert_close(*got, *want, n as f64, &format!("n = {n}, bin {k}"));
        }
    }
}

#[test]
fn inverse_transform_round_trips() {
    let mut rng = Rng::new(0x1234_5678_9ABC_DEF0);
    for log_n in 0..=10u32 {
        let n = 1usize << log_n;
        let input: Vec<Complex> = (0..n)
            .map(|_| Complex::new(rng.next_f64(), rng.next_f64()))
            .collect();
        let restored = ifft(&fft(&input).unwrap()).unwrap();
        for (index, (got, want)) in restored.iter().zip(input.iter()).enumerate() {
            assert_close(*got, *want, 1.0, &format!("n = {n}, sample {index}"));
        }
    }
}

#[test]
fn known_small_transforms_are_exact() {
    // [1, 2, 3, 4] -> [10, -2 + 2i, -2, -2 - 2i]. Derived by hand from the
    // definition; every term is a quarter-turn so the answer is exact.
    let input: Vec<Complex> = [1.0, 2.0, 3.0, 4.0]
        .iter()
        .map(|v| Complex::new(*v, 0.0))
        .collect();
    let out = fft(&input).unwrap();
    assert_eq!(out[0], Complex::new(10.0, 0.0));
    assert_eq!(out[1], Complex::new(-2.0, 2.0));
    assert_eq!(out[2], Complex::new(-2.0, 0.0));
    assert_eq!(out[3], Complex::new(-2.0, -2.0));

    // A unit impulse at the origin transforms to the constant 1.
    let mut impulse = vec![Complex::ZERO; 16];
    impulse[0] = Complex::new(1.0, 0.0);
    for bin in fft(&impulse).unwrap() {
        assert_eq!(bin, Complex::new(1.0, 0.0));
    }

    // A constant transforms to an impulse of height N at bin 0. The vanishing
    // bins are exact because every butterfly difference is exactly zero.
    let constant = vec![Complex::new(2.5, 0.0); 32];
    let out = fft(&constant).unwrap();
    assert_eq!(out[0], Complex::new(80.0, 0.0));
    for bin in &out[1..] {
        assert_eq!(bin.re, 0.0);
        assert_eq!(bin.im, 0.0);
    }
}

#[test]
fn real_transform_matches_the_complex_one() {
    let mut rng = Rng::new(0xDEAD_BEEF_CAFE_F00D);
    for log_n in 1..=10u32 {
        let n = 1usize << log_n;
        let values: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
        let packed = real_fft(&values).expect("power-of-two length");
        assert_eq!(packed.len(), n / 2 + 1);

        let as_complex: Vec<Complex> = values.iter().map(|v| Complex::new(*v, 0.0)).collect();
        let full = naive_dft(&as_complex);
        for (k, bin) in packed.iter().enumerate() {
            assert_close(*bin, full[k], n as f64, &format!("n = {n}, real bin {k}"));
        }
        // The two purely real bins are returned with an exact zero imaginary
        // part, not merely a small one.
        assert_eq!(packed[0].im, 0.0);
        assert_eq!(packed[n / 2].im, 0.0);

        let restored = real_ifft(&packed, n).unwrap();
        for (index, (got, want)) in restored.iter().zip(values.iter()).enumerate() {
            assert!(
                (got - want).abs() <= 1e-12,
                "n = {n}, sample {index}: {got} != {want}"
            );
        }
    }
}

#[test]
fn real_transform_end_bins_are_the_sums() {
    // X_0 is the plain sum and X_{N/2} the alternating sum, both by definition.
    let values: Vec<f64> = (0..8).map(|i| (i as f64) * 0.5 - 1.0).collect();
    let packed = real_fft(&values).unwrap();
    let sum: f64 = values.iter().sum();
    let alternating: f64 = values
        .iter()
        .enumerate()
        .map(|(i, v)| if i % 2 == 0 { *v } else { -*v })
        .sum();
    assert!((packed[0].re - sum).abs() <= 1e-12);
    assert!((packed[4].re - alternating).abs() <= 1e-12);
}

#[test]
fn single_point_transforms_are_the_identity() {
    // evergreen's `DIF<0>` specialisations do nothing, so a one-point transform
    // returns its input; the port reproduces that rather than erroring.
    let packed = real_fft(&[3.25]).unwrap();
    assert_eq!(packed, vec![Complex::new(3.25, 0.0)]);
    assert_eq!(real_ifft(&packed, 1).unwrap(), vec![3.25]);
    assert_eq!(
        fft(&[Complex::new(1.0, -2.0)]).unwrap()[0],
        Complex::new(1.0, -2.0)
    );
}

#[test]
fn the_packed_real_transform_still_requires_a_power_of_two() {
    // Its N/4 unpacking loop needs one, so this refusal is a property of the
    // algorithm and survives the move to rustfft.
    for bad_length in [3usize, 5, 6, 7, 100, 1000] {
        let values = vec![1.0; bad_length];
        assert!(matches!(real_fft(&values), Err(Error::InvalidValue(_))));
        assert!(matches!(
            real_ifft(&[Complex::ZERO], bad_length),
            Err(Error::InvalidValue(_))
        ));
    }
    assert!(matches!(fft(&[]), Err(Error::InvalidValue(_))));
    assert!(matches!(real_fft(&[]), Err(Error::InvalidValue(_))));
}

/// The complex transform used to refuse these lengths. That refusal guarded
/// against evergreen, which rounds `log2(len)` with its shape assertion compiled
/// out of release builds and so transforms a different number of points in
/// silence. rustfft transforms any length correctly, so the guard is gone and
/// the answer is checked against the naive DFT instead - including prime
/// lengths, which rustfft handles through Bluestein's algorithm.
#[test]
fn complex_transforms_of_any_length_match_the_naive_dft() {
    let mut rng = Rng::new(0x5eed_f00d);
    for n in [1usize, 2, 3, 5, 6, 7, 97, 100, 257, 1000] {
        let input: Vec<Complex> = (0..n)
            .map(|_| Complex::new(rng.next_f64(), rng.next_f64()))
            .collect();
        let expected = naive_dft(&input);
        let got = fft(&input).unwrap();
        for (k, (a, b)) in got.iter().zip(&expected).enumerate() {
            assert_close(*a, *b, n as f64, &format!("n={n} bin {k}"));
        }
        // And the inverse returns the input.
        let back = ifft(&got).unwrap();
        for (j, (a, b)) in back.iter().zip(&input).enumerate() {
            assert_close(*a, *b, 1.0, &format!("n={n} round trip {j}"));
        }
    }
}

#[test]
fn oversized_and_malformed_transforms_are_refused() {
    assert!(matches!(
        real_ifft(&[Complex::ZERO; 3], MAX_LEN * 2),
        Err(Error::InvalidRange(_))
    ));
    // An inverse real transform needs exactly N/2 + 1 bins.
    assert!(matches!(
        real_ifft(&[Complex::ZERO; 3], 8),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        real_fft(&[1.0, f64::NAN, 2.0, 3.0]),
        Err(Error::InvalidValue(_))
    ));
    assert!(matches!(
        real_ifft(
            &[
                Complex::new(f64::INFINITY, 0.0),
                Complex::ZERO,
                Complex::ZERO
            ],
            4
        ),
        Err(Error::InvalidValue(_))
    ));
}

#[test]
fn in_place_and_by_value_transforms_agree_bit_for_bit() {
    let mut rng = Rng::new(0x5DEE_CE66_D000_0001);
    let input: Vec<Complex> = (0..64)
        .map(|_| Complex::new(rng.next_f64(), rng.next_f64()))
        .collect();
    let by_value = fft(&input).unwrap();
    let mut in_place = input.clone();
    fft_in_place(&mut in_place).unwrap();
    assert_eq!(by_value, in_place);
}

#[test]
fn convolution_theorem_holds() {
    // An independent invariant: the inverse transform of a pointwise product is
    // the circular convolution. This is what the kernel density estimator
    // relies on, so it is worth asserting directly rather than inferring.
    let n = 16usize;
    let mut rng = Rng::new(0xABCD_0123_4567_89EF);
    let a: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();
    let b: Vec<f64> = (0..n).map(|_| rng.next_f64()).collect();

    let fa = real_fft(&a).unwrap();
    let fb = real_fft(&b).unwrap();
    let product: Vec<Complex> = fa.iter().zip(fb.iter()).map(|(x, y)| *x * *y).collect();
    let got = real_ifft(&product, n).unwrap();

    for (index, value) in got.iter().enumerate() {
        let mut expected = 0.0;
        for (j, aj) in a.iter().enumerate() {
            expected += aj * b[(index + n - j) % n];
        }
        assert!(
            (value - expected).abs() <= 1e-12,
            "position {index}: {value} != {expected}"
        );
    }
}
