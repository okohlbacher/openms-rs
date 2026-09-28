// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: OpenMS Rust contributors $
//
// No GNU C Library source was used: the instruction sequence and the constants
// below were read from the disassembly and `.rodata` of the reference build's
// `libm.so.6` (ibminode06, Ubuntu glibc 2.39-0ubuntu8.9, sha256
// fce00b6f25f459cf4ae0b7fae4257909a1a8a86f9149e7c029a5d617baf1ccd0), function
// `__logf_fma` (the ifunc choice on FMA hosts) at 0x7de50, table at 0xb7d40.
// The algorithm is Arm's optimized-routines `logf` (MIT OR Apache-2.0 WITH
// LLVM-exception), as in `crate::analysis::feature_finder_picked::glibc_libm`.

//! The single-precision natural logarithm of the reference build.
//!
//! `ElutionModelFitter::fitElutionModels` evaluates its imputation line at
//! `log((*it)->getIntensity())`, and `getIntensity()` returns a `float`, so
//! C++ overload resolution calls `std::log(float)`: the C library's `logf`.
//! A correctly rounded or host `logf` can differ from glibc's (which is not
//! correctly rounded) in the last bit, and that bit moves an imputed intensity
//! by several `float` ulps; this transcription makes the value identical to the
//! Linux x86_64 Release build on every host. `f64::mul_add` is the fused
//! multiply-add the reference instructions use.

/// `(invc, logc)` for the 16 subintervals, bit patterns from `.rodata`.
const TABLE: [(u64, u64); 16] = [
    (0x3ff6_61ec_79f8_f3be, 0xbfd5_7bf7_808c_aade),
    (0x3ff5_71ed_4aaf_883d, 0xbfd2_bef0_a7c0_6ddb),
    (0x3ff4_9539_f0f0_10b0, 0xbfd0_1eae_7f51_3a67),
    (0x3ff3_c995_b0b8_0385, 0xbfcb_31d8_a682_24e9),
    (0x3ff3_0d19_0c88_64a5, 0xbfc6_574f_0ac0_7758),
    (0x3ff2_5e22_7b0b_8ea0, 0xbfc1_aa2b_c79c_8100),
    (0x3ff1_bb4a_4a1a_343f, 0xbfba_4e76_ce8c_0e5e),
    (0x3ff1_2358_f08a_e5ba, 0xbfb1_973c_5a61_1ccc),
    (0x3ff0_953f_4199_00a7, 0xbfa2_52f4_38e1_0c1e),
    (0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000),
    (0x3fee_608c_fd9a_47ac, 0x3faa_a5aa_5df2_5984),
    (0x3fec_a4b3_1f02_6aa0, 0x3fbc_5e53_aa36_2eb4),
    (0x3feb_2036_576a_fce6, 0x3fc5_26e5_7720_db08),
    (0x3fe9_c2d1_63a1_aa2d, 0x3fcb_c286_0d22_4770),
    (0x3fe8_86e6_0378_41ed, 0x3fd1_058b_c8a0_7ee1),
    (0x3fe7_67dc_f553_4862, 0x3fd4_0430_57b6_ee09),
];
/// `Ln2` (0xb7e40).
const LN2: u64 = 0x3fe6_2e42_fefa_39ef;
/// Polynomial coefficients `A[0]` (0xb7e48), `A[1]` (0xb7e50), `A[2]` (0xb7e58).
const A0: u64 = 0xbfd0_0ea3_48b8_8334;
const A1: u64 = 0x3fd5_575b_0be0_0b6a;
const A2: u64 = 0xbfdf_fffe_f20a_4123;

/// `logf(x)` as the reference build's `__logf_fma` computes it.
///
/// Zero gives negative infinity, positive infinity itself, a negative or NaN
/// argument NaN (`(x - x) / (x - x)`, whose sign and payload are the host's),
/// and subnormals are scaled by 2^23 first, as the reference does.
pub(crate) fn logf(x: f32) -> f32 {
    let mut ix = x.to_bits();
    if ix == 0x3f80_0000 {
        return 0.0;
    }
    if ix.wrapping_sub(0x0080_0000) >= 0x7f80_0000 - 0x0080_0000 {
        if ix.wrapping_mul(2) == 0 {
            return f32::NEG_INFINITY;
        }
        if ix == 0x7f80_0000 {
            return x;
        }
        if ix & 0x8000_0000 != 0 || ix.wrapping_mul(2) >= 0xff00_0000 {
            return (x - x) / (x - x);
        }
        // Subnormal: normalise.
        ix = (x * f32::from_bits(0x4b00_0000)).to_bits();
        ix = ix.wrapping_sub(23 << 23);
    }
    let tmp = ix.wrapping_sub(0x3f33_0000);
    let i = ((tmp >> 19) & 15) as usize;
    let k = (tmp as i32) >> 23;
    let iz = ix.wrapping_sub(tmp & 0xff80_0000);
    let (invc, logc) = (f64::from_bits(TABLE[i].0), f64::from_bits(TABLE[i].1));
    let z = f64::from(f32::from_bits(iz));
    let y0 = f64::from(k).mul_add(f64::from_bits(LN2), logc);
    let r = z.mul_add(invc, -1.0);
    let y = r.mul_add(f64::from_bits(A1), f64::from_bits(A2));
    let r2 = r * r;
    let r_plus = r + y0;
    let y = r2.mul_add(f64::from_bits(A0), y);
    r2.mul_add(y, r_plus) as f32
}

#[cfg(test)]
mod tests {
    use super::logf;

    #[test]
    fn exact_points_and_special_values() {
        assert_eq!(logf(1.0), 0.0);
        assert_eq!(logf(0.0), f32::NEG_INFINITY);
        assert_eq!(logf(f32::INFINITY), f32::INFINITY);
        assert!(logf(-1.0).is_nan());
        assert!(logf(f32::NAN).is_nan());
        // Within one ulp of the correctly rounded value over a sweep.
        for k in 1..2000u32 {
            let x = f32::from_bits(0x0010_0000 + k * 0x0010_3f1d);
            let exact = f64::from(x).ln();
            let got = f64::from(logf(x));
            let ulp = f64::from(f32::from_bits((got as f32).to_bits() + 1)) - got;
            assert!((got - exact).abs() <= ulp.abs() + f64::EPSILON, "{x}");
        }
    }
}
