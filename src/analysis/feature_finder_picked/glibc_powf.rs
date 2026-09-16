// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause AND MIT
// $Maintainer: OpenMS Rust contributors $
//
// The algorithm and the two tables below are ported from Arm's
// optimized-routines (`math/powf.c`, `math/powf_log2_data.c`,
// `math/exp2f_data.c`; the same code as musl's `src/math/powf.c`), which carry
// this notice:
//
//   Copyright (c) 2017-2018, Arm Limited.
//   SPDX-License-Identifier: MIT OR Apache-2.0 WITH LLVM-exception
//
//   Permission is hereby granted, free of charge, to any person obtaining a
//   copy of this software and associated documentation files (the
//   "Software"), to deal in the Software without restriction, including
//   without limitation the rights to use, copy, modify, merge, publish,
//   distribute, sublicense, and/or sell copies of the Software, and to permit
//   persons to whom the Software is furnished to do so, subject to the
//   following conditions:
//
//   The above copyright notice and this permission notice shall be included
//   in all copies or substantial portions of the Software.
//
//   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
//   NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
//   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
//   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
//   USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// No GNU C Library source was used: the instruction sequence below was read
// from the disassembly of the reference build's `libm.so.6`.

//! The `powf` of the picked feature finder's reference build
//! (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`, the overall seed score).
//!
//! `FeatureFinderAlgorithmPicked::run_` computes the overall score of step 3.2
//! as `std::pow(float, float)` (`FeatureFinderAlgorithmPicked.cpp:506`), the
//! C library's `powf`. That function is not correctly rounded, and which
//! values it misrounds depends on the library and on the CPU. The Linux x86_64
//! Release build `openms4-release-bc9cc12-c19e494-174b576` binds
//! `powf@GLIBC_2.27` of the host's GNU C Library 2.39 (Ubuntu
//! `2.39-0ubuntu8.9`, `libm.so.6` sha256 `fce00b6f...`), an indirect function
//! that selects `__powf_fma` on a CPU with FMA and AVX2, such as the AMD EPYC
//! 7763 of the reference node (`../oracle/ffap-complete-fix1`: `powf_probe
//! where` resolves `powf` to offset `0x7df50`, which the installed debug
//! symbols name `__powf_fma`).
//!
//! That function is Arm's optimized-routines algorithm: `log2(x)` from a
//! 16-entry table and a degree-5 polynomial, `y * log2(x)`, and `2^t` from a
//! 32-entry table and a cubic, all in `double`, rounded once to `float`. The
//! FMA variant is the same C code compiled with `-mfma -mavx2`, and the
//! compiler fused every `a * b + c` into one `vfmadd` instruction (the
//! disassembly is `../oracle/ffap-complete-fix1/logs/powf_fma_disasm.txt`):
//! [`f64::mul_add`] is that instruction. The special cases follow the same
//! disassembly, including `__math_may_uflowf` for `y * log2(x)` in
//! `[-150, -149)` and the round-to-nearest outcome of the overflow rounding
//! check; the tables equal the bytes of `__powf_log2_data` and
//! `__exp2f_data` in that `libm.so.6`.
//!
//! NaN results carry the bits the SSE instructions give them: the first NaN
//! operand of the instruction, quieted, or the default NaN `0xffc00000` of an
//! invalid operation.
//!
//! Evidence (`../oracle/ffap-complete-fix1`, executed on ibminode06 through
//! `dlsym`): every binary32 `x` with `y = 1.0f / 3.0f`, the exponent the
//! algorithm uses, in 256 digests; `2^26` generated pairs in each of four sets
//! and the first 256 of each; and every pair of a 42-value special grid. The
//! unit tests replay the grid, the rows and `2^20`-pair digests.

/// `__powf_log2_data.tab`: `(invc, logc)` for each of the 16 subintervals.
const LOG2_TAB: [(u64, u64); 16] = [
    (0x3ff6_61ec_79f8_f3be, 0xbfde_fec6_5b96_3019),
    (0x3ff5_71ed_4aaf_883d, 0xbfdb_0b68_32d4_fca4),
    (0x3ff4_9539_f0f0_10b0, 0xbfd7_418b_0a1f_b77b),
    (0x3ff3_c995_b0b8_0385, 0xbfd3_9de9_1a6d_cf7b),
    (0x3ff3_0d19_0c88_64a5, 0xbfd0_1d9b_f3f2_b631),
    (0x3ff2_5e22_7b0b_8ea0, 0xbfc9_7c1d_1b3b_7af0),
    (0x3ff1_bb4a_4a1a_343f, 0xbfc2_f9e3_93af_3c9f),
    (0x3ff1_2358_f08a_e5ba, 0xbfb9_60cb_bf78_8d5c),
    (0x3ff0_953f_4199_00a7, 0xbfaa_6f9d_b647_5fce),
    (0x3ff0_0000_0000_0000, 0x0000_0000_0000_0000),
    (0x3fee_608c_fd9a_47ac, 0x3fb3_38ca_9f24_f53d),
    (0x3fec_a4b3_1f02_6aa0, 0x3fc4_76a9_5438_91ba),
    (0x3feb_2036_576a_fce6, 0x3fce_840b_4ac4_e4d2),
    (0x3fe9_c2d1_63a1_aa2d, 0x3fd4_0645_f0c6_651c),
    (0x3fe8_86e6_0378_41ed, 0x3fd8_8e9c_2c1b_9ff8),
    (0x3fe7_67dc_f553_4862, 0x3fdc_e0a4_4eb1_7bcc),
];

/// `__powf_log2_data.poly`: the coefficients `A[0]` to `A[4]`.
const LOG2_POLY: [u64; 5] = [
    0x3fd2_7616_c949_6e0b,
    0xbfd7_1969_a075_c67a,
    0x3fde_c70a_6ca7_badd,
    0xbfe7_1547_48be_f6c8,
    0x3ff7_1547_652a_b82b,
];

/// `__exp2f_data.tab`: `2^(i/32)` with `i << 47` subtracted, as bits.
const EXP2_TAB: [u64; 32] = [
    0x3ff0_0000_0000_0000,
    0x3fef_d9b0_d315_8574,
    0x3fef_b558_6cf9_890f,
    0x3fef_9301_d012_5b51,
    0x3fef_72b8_3c7d_517b,
    0x3fef_5487_3168_b9aa,
    0x3fef_387a_6e75_6238,
    0x3fef_1e9d_f51f_dee1,
    0x3fef_06fe_0a31_b715,
    0x3fee_f1a7_373a_a9cb,
    0x3fee_dea6_4c12_3422,
    0x3fee_ce08_6061_892d,
    0x3fee_bfda_d536_2a27,
    0x3fee_b42b_569d_4f82,
    0x3fee_ab07_dd48_5429,
    0x3fee_a47e_b03a_5585,
    0x3fee_a09e_667f_3bcd,
    0x3fee_9f75_e8ec_5f74,
    0x3fee_a114_73eb_0187,
    0x3fee_a589_994c_ce13,
    0x3fee_ace5_422a_a0db,
    0x3fee_b737_b0cd_c5e5,
    0x3fee_c491_82a3_f090,
    0x3fee_d503_b23e_255d,
    0x3fee_e89f_995a_d3ad,
    0x3fee_ff76_f2fb_5e47,
    0x3fef_199b_dd85_529c,
    0x3fef_3720_dcef_9069,
    0x3fef_5818_dcfb_a487,
    0x3fef_7c97_337b_9b5f,
    0x3fef_a4af_a2a4_90da,
    0x3fef_d076_5b6e_4540,
];

/// `__exp2f_data.shift_scaled`, `0x1.8p52 / 32`.
const EXP2_SHIFT: u64 = 0x42e8_0000_0000_0000;

/// `__exp2f_data.poly`: `C[0]` to `C[2]`.
const EXP2_POLY: [u64; 3] = [
    0x3fac_6af8_4b91_2394,
    0x3fce_bfce_50fa_c4f3,
    0x3fe6_2e42_ff0c_52d6,
];

/// Subinterval offset of the logarithm, `OFF`.
const OFF: u32 = 0x3f33_0000;

/// The sign of a negative result in the table index, `SIGN_BIAS`.
const SIGN_BIAS: u32 = 1 << (5 + 11);

/// `y * log2(x)` above which the result overflows: `0x1.fffffffd1d571p+6`.
const OVERFLOW_BOUND: u64 = 0x405f_ffff_ffd1_d571;

/// `y * log2(x)` at or below which the result underflows to zero: `-150`.
const UNDERFLOW_BOUND: u64 = 0xc062_c000_0000_0000;

/// `y * log2(x)` below which the result may underflow: `-149`.
const MAY_UNDERFLOW_BOUND: u64 = 0xc062_a000_0000_0000;

/// SSE's default ("real indefinite") `float` NaN.
const DEFAULT_NAN: u32 = 0xffc0_0000;

fn quiet(x: f32) -> f32 {
    f32::from_bits(x.to_bits() | 0x0040_0000)
}

/// The NaN rule of a two-operand SSE `float` instruction whose first source
/// operand is `first`: that operand's NaN, quieted; else the second's; else
/// the default NaN for an invalid operation.
fn sse(first: f32, second: f32, result: f32) -> f32 {
    if !result.is_nan() {
        result
    } else if first.is_nan() {
        quiet(first)
    } else if second.is_nan() {
        quiet(second)
    } else {
        f32::from_bits(DEFAULT_NAN)
    }
}

/// `vaddss`: `first + second`.
fn add(first: f32, second: f32) -> f32 {
    sse(first, second, first + second)
}

/// `mulss`: `first * second`, with SSE's NaN rule.
pub(crate) fn mul(first: f32, second: f32) -> f32 {
    sse(first, second, first * second)
}

/// `vdivss`: `first / second`.
fn div(first: f32, second: f32) -> f32 {
    sse(first, second, first / second)
}

/// `zeroinfnan`: `x` is zero, infinite or NaN.
fn zero_inf_nan(bits: u32) -> bool {
    bits.wrapping_mul(2).wrapping_sub(1) >= 0xfeff_ffff
}

/// `issignalingf_inline`.
fn is_signaling(x: f32) -> bool {
    ((x.to_bits() ^ 0x0040_0000) & 0x7fff_ffff) > 0x7fc0_0000
}

/// `checkint`: 0 when `y` is not an integer, 1 when it is odd, 2 when even.
fn check_int(iy: u32) -> u32 {
    let exponent = (iy >> 23) & 0xff;
    if exponent < 0x7f {
        return 0;
    }
    if exponent > 0x7f + 23 {
        return 2;
    }
    let unit = 1u32 << (0x7f + 23 - exponent);
    if iy & (unit - 1) != 0 {
        0
    } else if iy & unit != 0 {
        1
    } else {
        2
    }
}

/// `xflowf`: `(sign ? -y : y) * y`, rounded.
fn xflow(sign: u32, y: f32) -> f32 {
    (if sign != 0 { -y } else { y }) * y
}

/// `__math_invalidf` for a finite `x`: `(x - x) / (x - x)`, where `x - x` is
/// `+0.0`, so the default NaN.
fn invalid(_x: f32) -> f32 {
    div(0.0, 0.0)
}

/// `__math_divzerof`: `±1 / 0`.
fn divide_by_zero(sign: u32) -> f32 {
    if sign != 0 {
        f32::NEG_INFINITY
    } else {
        f32::INFINITY
    }
}

/// `log2_inline`: `log2(x)` for the biased bits `ix` of a positive normalised
/// `x`, with the fused multiply-adds of `__powf_fma`.
fn log2(ix: u32) -> f64 {
    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> 19) & 0xf) as usize;
    let top = tmp & 0xff80_0000;
    let iz = ix.wrapping_sub(top);
    let k = (top as i32) >> 23;
    let (invc, logc) = LOG2_TAB[i];
    let z = f64::from(f32::from_bits(iz));
    let a = LOG2_POLY.map(f64::from_bits);
    let r = z.mul_add(f64::from_bits(invc), -1.0);
    let y0 = f64::from(k) + f64::from_bits(logc);
    let high = r.mul_add(a[0], a[1]);
    let p = r.mul_add(a[2], a[3]);
    let r2 = r * r;
    let q = r.mul_add(a[4], y0);
    let r4 = r2 * r2;
    let q = r2.mul_add(p, q);
    high.mul_add(r4, q)
}

/// `exp2_inline`: `2^xd`, negated when `sign_bias` says so, rounded to
/// `float`, with the fused multiply-adds of `__powf_fma`.
fn exp2(xd: f64, sign_bias: u32) -> f32 {
    let shift = f64::from_bits(EXP2_SHIFT);
    let c = EXP2_POLY.map(f64::from_bits);
    let kd = xd + shift;
    let ki = kd.to_bits();
    let kd = kd - shift;
    let r = xd - kd;
    let t =
        EXP2_TAB[(ki & 0x1f) as usize].wrapping_add(ki.wrapping_add(u64::from(sign_bias)) << 47);
    let s = f64::from_bits(t);
    let z = r.mul_add(c[0], c[1]);
    let r2 = r * r;
    let y = r.mul_add(c[2], 1.0);
    let y = z.mul_add(r2, y);
    (y * s) as f32
}

/// `x` to the power `y` as `__powf_fma` of the reference build's GNU C
/// Library 2.39 computes it, NaN bits included (module documentation).
pub(crate) fn powf(x: f32, y: f32) -> f32 {
    let mut ix = x.to_bits();
    let iy = y.to_bits();
    let mut sign_bias = 0u32;
    if ix.wrapping_sub(0x0080_0000) >= 0x7f80_0000 - 0x0080_0000 || zero_inf_nan(iy) {
        if zero_inf_nan(iy) {
            if iy.wrapping_mul(2) == 0 {
                return if is_signaling(x) { add(x, y) } else { 1.0 };
            }
            if ix == 0x3f80_0000 {
                return if is_signaling(y) { add(x, y) } else { 1.0 };
            }
            let (x2, y2) = (ix.wrapping_mul(2), iy.wrapping_mul(2));
            if x2 > 0xff00_0000 || y2 > 0xff00_0000 {
                return add(x, y);
            }
            if x2 == 0x7f00_0000 {
                return 1.0;
            }
            if (x2 < 0x7f00_0000) == (iy & 0x8000_0000 == 0) {
                return 0.0;
            }
            return mul(y, y);
        }
        if zero_inf_nan(ix) {
            let mut x2 = mul(x, x);
            if ix & 0x8000_0000 != 0 && check_int(iy) == 1 {
                x2 = -x2;
                sign_bias = 1;
            }
            if ix.wrapping_mul(2) == 0 && iy & 0x8000_0000 != 0 {
                return divide_by_zero(sign_bias);
            }
            return if iy & 0x8000_0000 != 0 {
                div(1.0, x2)
            } else {
                x2
            };
        }
        // `x` and `y` are finite and non-zero.
        if ix & 0x8000_0000 != 0 {
            match check_int(iy) {
                0 => return invalid(x),
                1 => sign_bias = SIGN_BIAS,
                _ => {}
            }
            ix &= 0x7fff_ffff;
        }
        if ix < 0x0080_0000 {
            // A subnormal `x`, normalised so that its exponent is negative.
            ix = (x * f32::from_bits(0x4b00_0000)).to_bits() & 0x7fff_ffff;
            ix = ix.wrapping_sub(23 << 23);
        }
    }
    let ylogx = f64::from(y) * log2(ix);
    if ((ylogx.to_bits() >> 47) & 0xffff) >= (126.0_f64.to_bits() >> 47) {
        if ylogx > f64::from_bits(OVERFLOW_BOUND) {
            return xflow(sign_bias, f32::from_bits(0x7000_0000));
        }
        // Between `0x1.fffffffa3aae2p+6` and the overflow bound the library
        // checks the rounding mode (`1.0f + 0x1p-25f != 1.0f`), which is false
        // in the round-to-nearest mode every Rust program runs in, and goes on
        // to `exp2_inline`, as below.
        if ylogx <= f64::from_bits(UNDERFLOW_BOUND) {
            return xflow(sign_bias, f32::from_bits(0x1000_0000));
        }
        if ylogx < f64::from_bits(MAY_UNDERFLOW_BOUND) {
            return xflow(sign_bias, f32::from_bits(0x1a20_0000));
        }
    }
    exp2(ylogx, sign_bias)
}

#[cfg(test)]
mod tests {
    use super::*;

    const THIRD: u32 = 0x3eaa_aaab;

    fn fixture() -> String {
        std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/feature_finder_picked/glibc_powf_probe.tsv"
        ))
        .unwrap()
    }

    fn hex(text: &str) -> u32 {
        u32::from_str_radix(text, 16).unwrap()
    }

    fn splitmix64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    const EXPONENTS: [u32; 16] = [
        0x3eaa_aaab,
        0x3f00_0000,
        0x4000_0000,
        0x4040_0000,
        0xbf80_0000,
        0xc000_0000,
        0x3a83_126f,
        0x44fa_0000,
        0x3dcc_cccd,
        0xbeaa_aaab,
        0x42c8_0000,
        0xc2c8_0000,
        0x3f80_0000,
        0x4120_0000,
        0x3fc0_0000,
        0x3eaa_aaaa,
    ];

    /// The driver's generator `draw` (`powf_probe.c`).
    fn draw(set: u64, state: &mut u64) -> (u32, u32) {
        let d = splitmix64(state);
        match set {
            1 => (d as u32, (d >> 32) as u32),
            2 => (d as u32, EXPONENTS[(d >> 60) as usize]),
            3 => {
                let ex = (((d >> 32) & 0xff) % 81) as u32 + 127 - 40;
                let ey = (((d >> 40) & 0xff) % 21) as u32 + 127 - 10;
                (
                    (ex << 23) | (d as u32 & 0x7f_ffff),
                    (d as u32 & 0x8000_0000) | (ey << 23) | (((d >> 48) as u32) << 7),
                )
            }
            _ => (d as u32 & 0x7fff_ffff, THIRD),
        }
    }

    fn fnv32(mut h: u64, v: u32) -> u64 {
        for byte in v.to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        h
    }

    #[test]
    fn the_resolved_function_is_the_fma_variant_of_glibc_2_39() {
        let text = fixture();
        let facts: Vec<&str> = text.lines().filter(|l| l.starts_with("where\t")).collect();
        assert_eq!(
            facts,
            [
                "where\tglibc\t2.39",
                "where\tobject\t/lib/x86_64-linux-gnu/libm.so.6",
                "where\toffset\t0x7df50",
            ]
        );
    }

    #[test]
    fn every_special_pair_matches_the_executed_library() {
        let text = fixture();
        let mut count = 0;
        for line in text.lines().filter(|l| l.starts_with("grid\t")) {
            let cells: Vec<&str> = line.split('\t').collect();
            let (x, y, r) = (hex(cells[1]), hex(cells[2]), hex(cells[3]));
            let got = powf(f32::from_bits(x), f32::from_bits(y)).to_bits();
            assert_eq!(got, r, "powf({x:08x}, {y:08x})");
            count += 1;
        }
        assert_eq!(count, 42 * 42);
    }

    #[test]
    fn generated_pairs_match_the_executed_library() {
        let text = fixture();
        let pairs: Vec<Vec<&str>> = text
            .lines()
            .filter(|l| l.starts_with("pair\t"))
            .map(|l| l.split('\t').collect())
            .collect();
        assert_eq!(pairs.len(), 4 * 256);
        for cells in &pairs {
            let (x, y, r) = (hex(cells[2]), hex(cells[3]), hex(cells[4]));
            let got = powf(f32::from_bits(x), f32::from_bits(y)).to_bits();
            assert_eq!(got, r, "set {} powf({x:08x}, {y:08x})", cells[1]);
        }
        let digests: Vec<Vec<&str>> = text
            .lines()
            .filter(|l| l.starts_with("digest\t"))
            .map(|l| l.split('\t').collect())
            .collect();
        assert_eq!(digests.len(), 4);
        for (set, digest) in (1u64..).zip(&digests) {
            let n: u64 = digest[2].parse().unwrap();
            assert_eq!(n, 1 << 20);
            let mut state = set;
            let mut h = 0xcbf2_9ce4_8422_2325u64;
            for i in 0..n {
                let (x, y) = draw(set, &mut state);
                if i < 256 {
                    // The listed rows are the first draws of the same generator.
                    let row = &pairs[((set - 1) * 256 + i) as usize];
                    assert_eq!((hex(row[2]), hex(row[3])), (x, y), "set {}", digest[1]);
                }
                let r = powf(f32::from_bits(x), f32::from_bits(y)).to_bits();
                h = fnv32(fnv32(fnv32(h, x), y), r);
            }
            assert_eq!(format!("{h:016x}"), digest[3], "set {}", digest[1]);
        }
    }
}
