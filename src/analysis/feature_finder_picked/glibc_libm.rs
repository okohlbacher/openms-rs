// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause AND MIT
// $Maintainer: OpenMS Rust contributors $
//
// The algorithms of `exp` and `log` and their tables below are ported from
// Arm's optimized-routines (`math/exp.c`, `math/exp_data.c`, `math/log.c`,
// `math/log_data.c`, master 5288c42d; the same code as musl's
// `src/math/exp.c` and `src/math/log.c`), which carry this notice:
//
//   Copyright (c) 2018-2025, Arm Limited.
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
// No GNU C Library source was used: the instruction sequences below were read
// from the disassembly of the reference build's `libm.so.6`.

//! The `exp`, `log`, `atan` and `sqrt` of the picked feature finder's
//! reference build (`FEATUREFINDER/FeatureFinderAlgorithmPicked.h`: the
//! Gaussian and EGH trace fitters, their gnuplot and area formulas).
//!
//! `GaussTraceFitter` and `EGHTraceFitter` call the C library's `exp`, `log`,
//! `atan` and `sqrt` (`libOpenMS.so` calls `exp@plt`, `log@plt`, `atan@plt`
//! and `sqrt@plt` from both fitters; `../oracle/ffap-complete-fix3`,
//! `logs/trace_fitters_libm_calls.txt`). None of the transcendental ones is
//! correctly rounded, and which values they misround depends on the library
//! and, through indirect functions, on the CPU. The Linux x86_64 Release build
//! `openms4-release-bc9cc12-c19e494-174b576` binds `exp@GLIBC_2.29`,
//! `log@GLIBC_2.29` and `atan@GLIBC_2.2.5` of the host's GNU C Library 2.39
//! (Ubuntu `2.39-0ubuntu8.9`, `libm.so.6` sha256 `fce00b6f...`). On the AMD
//! EPYC 7763 of the reference node, which has FMA and AVX2, the exported
//! `exp` (`__GI___exp`) calls the indirect function `__ieee754_exp`, which
//! selects `__ieee754_exp_fma`; the exported `log` (`__log`) jumps to
//! `__ieee754_log`, which selects `__ieee754_log_fma`; `atan` selects
//! `__atan_fma` (the resolvers and the installed debug symbols,
//! `logs/libm_disasm.txt` and `logs/libm_disasm2.txt`). Both exported
//! wrappers only set `errno` and return the value unchanged.
//!
//! # `exp` and `log`: ported
//!
//! `__ieee754_exp_fma` and `__ieee754_log_fma` are Arm's optimized-routines
//! algorithms (128-entry tables; a degree-5 polynomial for `exp`, degree 6
//! for `log` and degree 12 near 1), compiled with `-mfma -mavx2`; the
//! compiler fused multiply-adds into `vfmadd`/`vfnmadd` instructions, and
//! [`f64::mul_add`] is that instruction. Which products are fused was read
//! from the disassembly and is noted at each step: `exp` fuses the rounding
//! `x * N/ln2 + shift` and the reduction, the polynomial, and the final
//! `scale + scale * tmp` except in the subnormal branch; `log` computes
//! `z / c - 1` as `fma(z, 1/c, -1)` (the `__FP_FAST_FMA` branch, which needs
//! no second table) and fuses the near-1 splitting `r + r * 2^27` and its
//! subtraction, so the split is not the one a compiler without contraction
//! would compute. The tables and constants equal the bytes of `__exp_data`,
//! `__log_data` and the literal pool in that `libm.so.6`
//! (`extract/check_tables.py`). The special cases follow the same
//! disassembly: `__math_oflow`/`__math_uflow` return `0x1p769 * 0x1p769` and
//! `0x1p-767 * 0x1p-767`, and `__math_invalid` computes `(x - x) / (x - x)`.
//!
//! NaN results carry the bits the SSE instructions give them: the first NaN
//! operand of the instruction, quieted, or the default NaN
//! `0xfff8000000000000` of an invalid operation (`log` of a negative number
//! or of `-inf`). No fused instruction ever sees a NaN.
//!
//! Evidence (`../oracle/ffap-complete-fix3`, `libm_probe`, executed on
//! ibminode06 through `dlsym`, twice, identical): an 80-value special grid
//! through each function; for each of thirteen generated sets (five for
//! `exp`, five for `log`, three for `atan`) of `2^26` inputs, 64 digests of
//! `2^20` inputs and the first 256 rows. The port equals the executed library
//! on every `exp` and `log` input (`port-harness`, Linux x86_64 and macOS
//! arm64); the unit tests replay the grid, the rows and the first digest of
//! every set.
//!
//! # `atan`: the host's
//!
//! `__atan_fma` is the IBM Accurate Mathematical Library's `atan` without
//! its multi-precision slow path. That library is licensed only under the
//! LGPL inside the GNU C Library; there is no MIT, Apache or fdlibm-style
//! upstream of the same algorithm (Arm optimized-routines has no scalar
//! double `atan`, and musl, fdlibm and CORE-MATH implement other algorithms).
//! Lead decision D10's fallback therefore applies: on x86_64 Linux with the
//! GNU C Library, [`atan`] calls the host's `atan` ([`f64::atan`]). That is
//! the reference build's function, and the result exact, only where the
//! host's library selects `__atan_fma`: GNU C Library 2.39 on a CPU with FMA,
//! as on the reference node (the port harness replays every `atan` input of
//! the probe on kim). Another GNU C Library version, or a CPU without FMA,
//! whose resolver selects another variant, is not measured. Everywhere else
//! (macOS, Windows, other architectures, other C libraries) it calls the
//! `libm` crate's `atan` (FreeBSD's `s_atan.c`), which gives the same value on
//! every such host; a NaN argument returns `x + x`, as `__atan_fma` does.
//!
//! How far that fallback departs from `__atan_fma`, measured by the round-3
//! verifier (`../oracle/ffc-numerics-v3`, `results/harness/node.txt`, `2^28`
//! inputs per set): it differs on 16,584,995 inputs `x` in `[0, 10]` (6.2 %),
//! on 4,174,993 inputs with `|x|` in `[2^-14, 2^15)` (1.6 %) and on 59,108
//! random bit patterns (0.02 %). Neither is correctly rounded: on 20,000
//! inputs of each of the first two ranges the reference misrounds 13 and 7,
//! the `libm` crate 1,219 and 324 (80-digit decimal check,
//! `extract/atan_rounding.py`), so a correctly rounded `atan` would depart
//! from the reference about 90 times less often, but it is still not the
//! reference algorithm and would need a new dependency. Only
//! `EGHTraceFitter::getArea` calls `atan`. The feature intensity is that area
//! divided by the window maximum and narrowed to `float`
//! (`FeatureFinderAlgorithmPicked.cpp:790`), which hides nearly every
//! last-bit difference: every fixture of this crate, the 67 returned EGH runs
//! with 767 features of `extended_stage.tsv.gz` included, is bit for bit on
//! macOS arm64 (see
//! `docs/EGH_TRACE_FITTER_SUPPORT.md`). The `double` area itself still
//! differs for a few percent of fits on such hosts (the round-3 verifier's
//! measurement), so the bound is a measured maximum over the listed
//! fixtures, not a guarantee.
//!
//! # `sqrt`
//!
//! `sqrt` is correctly rounded everywhere; [`sqrt`] only gives a NaN result
//! the bits `sqrtsd` gives it (`x86_64::sqrt` of the scoring module).

use crate::analysis::feature_finder_picked::scoring::x86_64;

/// `__exp_data.tab`: for each `k` in `0..128`, `T[k]` as bits and then the
/// bits of `H[k]` minus `(k << 52) / 128`, where `2^(k/128) ~= H[k] * (1 + T[k])`.
const EXP_TAB: [u64; 256] = [
    0x0000_0000_0000_0000,
    0x3ff0_0000_0000_0000,
    0x3c9b_3b4f_1a88_bf6e,
    0x3fef_f63d_a9fb_3335,
    0xbc71_6013_9cd8_dc5d,
    0x3fef_ec9a_3e77_8061,
    0xbc90_5e7a_1087_66d1,
    0x3fef_e315_e86e_7f85,
    0x3c8c_d252_3567_f613,
    0x3fef_d9b0_d315_8574,
    0xbc8b_ce80_23f9_8efa,
    0x3fef_d06b_29dd_f6de,
    0x3c60_f74e_61e6_c861,
    0x3fef_c745_1875_9bc8,
    0x3c90_a3e4_5b33_d399,
    0x3fef_be3e_cac6_f383,
    0x3c97_9aa6_5d83_7b6d,
    0x3fef_b558_6cf9_890f,
    0x3c8e_b51a_92fd_effc,
    0x3fef_ac92_2b72_47f7,
    0x3c3e_be3d_702f_9cd1,
    0x3fef_a3ec_32d3_d1a2,
    0xbc6a_0334_8990_6e0b,
    0x3fef_9b66_affe_d31b,
    0xbc95_5652_2a2f_bd0e,
    0x3fef_9301_d012_5b51,
    0xbc50_80ef_8c4e_ea55,
    0x3fef_8abd_c06c_31cc,
    0xbc91_c923_b9d5_f416,
    0x3fef_829a_aea9_2de0,
    0x3c80_d3e3_e95c_55af,
    0x3fef_7a98_c8a5_8e51,
    0xbc80_1b15_eaa5_9348,
    0x3fef_72b8_3c7d_517b,
    0xbc8f_1ff0_55de_323d,
    0x3fef_6af9_388c_8dea,
    0x3c8b_898c_3f13_53bf,
    0x3fef_635b_eb6f_cb75,
    0xbc96_d99c_7611_eb26,
    0x3fef_5be0_8404_5cd4,
    0x3c9a_ecf7_3e3a_2f60,
    0x3fef_5487_3168_b9aa,
    0xbc8f_e782_cb86_389d,
    0x3fef_4d50_22fc_d91d,
    0x3c8a_6f41_44a6_c38d,
    0x3fef_463b_8862_8cd6,
    0x3c80_7a05_b0e4_047d,
    0x3fef_3f49_917d_dc96,
    0x3c96_8efd_e3a8_a894,
    0x3fef_387a_6e75_6238,
    0x3c87_5e18_f274_487d,
    0x3fef_31ce_4fb2_a63f,
    0x3c80_472b_981f_e7f2,
    0x3fef_2b45_65e2_7cdd,
    0xbc96_b87b_3f71_085e,
    0x3fef_24df_e1f5_6381,
    0x3c82_f7e1_6d09_ab31,
    0x3fef_1e9d_f51f_dee1,
    0xbc3d_219b_1a6f_bffa,
    0x3fef_187f_d0da_d990,
    0x3c8b_3782_720c_0ab4,
    0x3fef_1285_a6e4_030b,
    0x3c6e_1492_89ce_cb8f,
    0x3fef_0caf_a93e_2f56,
    0x3c83_4d75_4db0_abb6,
    0x3fef_06fe_0a31_b715,
    0x3c86_4201_e2ac_744c,
    0x3fef_0170_fc4c_d831,
    0x3c8f_dd39_5dd3_f84a,
    0x3fee_fc08_b264_16ff,
    0xbc86_a380_3b8e_5b04,
    0x3fee_f6c5_5f92_9ff1,
    0xbc92_4aed_cc4b_5068,
    0x3fee_f1a7_373a_a9cb,
    0xbc99_07f8_1b51_2d8e,
    0x3fee_ecae_6d05_d866,
    0xbc71_d1e8_3e94_36d2,
    0x3fee_e7db_34e5_9ff7,
    0xbc99_1919_b3ce_1b15,
    0x3fee_e32d_c313_a8e5,
    0x3c85_9f48_a72a_4c6d,
    0x3fee_dea6_4c12_3422,
    0xbc93_1260_7a28_698a,
    0x3fee_da45_04ac_801c,
    0xbc58_a78f_4817_895b,
    0x3fee_d60a_21f7_2e2a,
    0xbc7c_2c9b_6749_9a1b,
    0x3fee_d1f5_d950_a897,
    0x3c43_63ed_60c2_ac11,
    0x3fee_ce08_6061_892d,
    0x3c96_6609_3b06_64ef,
    0x3fee_ca41_ed1d_0057,
    0x3c6e_cce1_daa1_0379,
    0x3fee_c6a2_b5c1_3cd0,
    0x3c93_ff8e_3f0f_1230,
    0x3fee_c32a_f0d7_d3de,
    0x3c76_90ce_bb7a_afb0,
    0x3fee_bfda_d536_2a27,
    0x3c93_1dbd_eb54_e077,
    0x3fee_bcb2_99fd_dd0d,
    0xbc8f_9434_0071_a38e,
    0x3fee_b9b2_769d_2ca7,
    0xbc87_decc_dc93_a349,
    0x3fee_b6da_a2cf_6642,
    0xbc78_dec6_bd0f_385f,
    0x3fee_b42b_569d_4f82,
    0xbc86_1246_ec7b_5cf6,
    0x3fee_b1a4_ca5d_920f,
    0x3c93_3505_18fd_d78e,
    0x3fee_af47_36b5_27da,
    0x3c7b_98b7_2f8a_9b05,
    0x3fee_ad12_d497_c7fd,
    0x3c90_63e1_e21c_5409,
    0x3fee_ab07_dd48_5429,
    0x3c34_c785_5019_c6ea,
    0x3fee_a926_8a59_46b7,
    0x3c94_32e6_2b64_c035,
    0x3fee_a76f_15ad_2148,
    0xbc8c_e44a_6199_769f,
    0x3fee_a5e1_b976_dc09,
    0xbc8c_33c5_3bef_4da8,
    0x3fee_a47e_b03a_5585,
    0xbc84_5378_892b_e9ae,
    0x3fee_a346_34cc_c320,
    0xbc93_cedd_7856_5858,
    0x3fee_a238_8255_2225,
    0x3c57_10aa_807e_1964,
    0x3fee_a155_d44c_a973,
    0xbc93_b3ef_bf5e_2228,
    0x3fee_a09e_667f_3bcd,
    0xbc6a_12ad_8734_b982,
    0x3fee_a012_750b_dabf,
    0xbc63_67ef_b86d_a9ee,
    0x3fee_9fb2_3c65_1a2f,
    0xbc80_dc3d_54e0_8851,
    0x3fee_9f7d_f951_9484,
    0xbc78_1f64_7e5a_3ecf,
    0x3fee_9f75_e8ec_5f74,
    0xbc86_ee4a_c08b_7db0,
    0x3fee_9f9a_48a5_8174,
    0xbc86_1932_1e55_e68a,
    0x3fee_9feb_5642_67c9,
    0x3c90_9ccb_5e09_d4d3,
    0x3fee_a069_4fde_5d3f,
    0xbc7b_32dc_b94d_a51d,
    0x3fee_a114_73eb_0187,
    0x3c94_ecfd_5467_c06b,
    0x3fee_a1ed_0130_c132,
    0x3c65_ebe1_abd6_6c55,
    0x3fee_a2f3_36cf_4e62,
    0xbc88_a1c5_2fb3_cf42,
    0x3fee_a427_543e_1a12,
    0xbc93_69b6_f13b_3734,
    0x3fee_a589_994c_ce13,
    0xbc80_5e84_3a19_ff1e,
    0x3fee_a71a_4623_c7ad,
    0xbc94_d450_d872_576e,
    0x3fee_a8d9_9b44_92ed,
    0x3c90_ad67_5b0e_8a00,
    0x3fee_aac7_d98a_6699,
    0x3c8d_b72f_c1f0_eab4,
    0x3fee_ace5_422a_a0db,
    0xbc65_b660_9cc5_e7ff,
    0x3fee_af32_16b5_448c,
    0x3c7b_f683_59f3_5f44,
    0x3fee_b1ae_9915_7736,
    0xbc93_091f_a71e_3d83,
    0x3fee_b45b_0b91_ffc6,
    0xbc5d_a9b8_8b6c_1e29,
    0x3fee_b737_b0cd_c5e5,
    0xbc6c_23f9_7c90_b959,
    0x3fee_ba44_cbc8_520f,
    0xbc92_4343_22f4_f9aa,
    0x3fee_bd82_9fde_4e50,
    0xbc85_ca6c_d766_8e4b,
    0x3fee_c0f1_70ca_07ba,
    0x3c71_affc_2b91_ce27,
    0x3fee_c491_82a3_f090,
    0x3c6d_d235_e10a_73bb,
    0x3fee_c863_19e3_2323,
    0xbc87_c504_2262_2263,
    0x3fee_cc66_7b5d_e565,
    0x3c8b_1c86_e3e2_31d5,
    0x3fee_d09b_ec4a_2d33,
    0xbc91_bbd1_d3bc_bb15,
    0x3fee_d503_b23e_255d,
    0x3c90_cc31_9cee_31d2,
    0x3fee_d99e_1330_b358,
    0x3c84_6984_6e73_5ab3,
    0x3fee_de6b_5579_fdbf,
    0xbc82_dfcd_978e_9db4,
    0x3fee_e36b_bfd3_f37a,
    0x3c8c_1a77_92cb_3387,
    0x3fee_e89f_995a_d3ad,
    0xbc90_7b8f_4ad1_d9fa,
    0x3fee_ee07_298d_b666,
    0xbc55_c3d9_56dc_aeba,
    0x3fee_f3a2_b84f_15fb,
    0xbc90_a40e_3da6_f640,
    0x3fee_f972_8de5_593a,
    0xbc68_d6f4_38ad_9334,
    0x3fee_ff76_f2fb_5e47,
    0xbc91_eee2_6b58_8a35,
    0x3fef_05b0_30a1_064a,
    0x3c74_ffd7_0a5f_ddcd,
    0x3fef_0c1e_904b_c1d2,
    0xbc91_bdfb_fa92_98ac,
    0x3fef_12c2_5bd7_1e09,
    0x3c73_6eae_30af_0cb3,
    0x3fef_199b_dd85_529c,
    0x3c8e_e332_5c9f_fd94,
    0x3fef_20ab_5fff_d07a,
    0x3c84_e08f_d109_59ac,
    0x3fef_27f1_2e57_d14b,
    0x3c63_cdaf_384e_1a67,
    0x3fef_2f6d_9406_e7b5,
    0x3c67_6b2c_6c92_1968,
    0x3fef_3720_dcef_9069,
    0xbc80_8a18_83cc_b5d2,
    0x3fef_3f0b_555d_c3fa,
    0xbc8f_ad5d_3fff_fa6f,
    0x3fef_472d_4a07_897c,
    0xbc90_0dae_3875_a949,
    0x3fef_4f87_080d_89f2,
    0x3c74_a385_a63d_07a7,
    0x3fef_5818_dcfb_a487,
    0xbc82_919e_2040_220f,
    0x3fef_60e3_16c9_8398,
    0x3c8e_5a50_d5c1_92ac,
    0x3fef_69e6_03db_3285,
    0x3c84_3a59_ac01_6b4b,
    0x3fef_7321_f301_b460,
    0xbc82_d521_07b4_3e1f,
    0x3fef_7c97_337b_9b5f,
    0xbc89_2ab9_3b47_0dc9,
    0x3fef_8646_14f5_a129,
    0x3c74_b604_603a_88d3,
    0x3fef_902e_e78b_3ff6,
    0x3c83_c5ec_519d_7271,
    0x3fef_9a51_fbc7_4c83,
    0xbc8f_f712_8fd3_91f0,
    0x3fef_a4af_a2a4_90da,
    0xbc8d_ae98_e223_747d,
    0x3fef_af48_2d8e_67f1,
    0x3c8e_c3bc_41aa_2008,
    0x3fef_ba1b_ee61_5a27,
    0x3c84_2b94_c3a9_eb32,
    0x3fef_c52b_376b_ba97,
    0x3c8a_64a9_31d1_85ee,
    0x3fef_d076_5b6e_4540,
    0xbc8e_37ba_e43b_e3ed,
    0x3fef_dbfd_ad9c_be14,
    0x3c77_893b_4d91_cd9d,
    0x3fef_e7c1_819e_90d8,
    0x3c53_05c1_4160_cc89,
    0x3fef_f3c2_2b8f_71f1,
];

/// `__log_data.tab`: `(invc, logc)` bits for each of the 128 subintervals.
const LOG_TAB: [(u64, u64); 128] = [
    (0x3ff7_34f0_c3e0_de9f, 0xbfd7_cc7f_79e6_9000),
    (0x3ff7_1378_6a2c_e91f, 0xbfd7_6fee_c20d_0000),
    (0x3ff6_f260_08fa_b5a0, 0xbfd7_13e3_1351_e000),
    (0x3ff6_d1a6_1f13_8c7d, 0xbfd6_b85b_3828_7800),
    (0x3ff6_b149_0bc5_b4d1, 0xbfd6_5d55_9080_7800),
    (0x3ff6_9147_332f_0cba, 0xbfd6_02d0_7618_0000),
    (0x3ff6_719f_1822_4223, 0xbfd5_a8ca_8690_9000),
    (0x3ff6_524f_99a5_1ed9, 0xbfd5_4f43_5603_5000),
    (0x3ff6_3356_aa8f_24c4, 0xbfd4_f637_c36b_4000),
    (0x3ff6_14b3_6b9d_dc14, 0xbfd4_9da7_fda8_5000),
    (0x3ff5_f664_52c6_5c4c, 0xbfd4_4592_3989_a800),
    (0x3ff5_d867_b591_2c4f, 0xbfd3_edf4_39b0_b800),
    (0x3ff5_babc_cb5b_90de, 0xbfd3_96ce_448f_7000),
    (0x3ff5_9d61_f2d9_1a78, 0xbfd3_401e_17bd_a000),
    (0x3ff5_8056_1246_5687, 0xbfd2_e9e2_ef46_8000),
    (0x3ff5_6397_cee7_6bd3, 0xbfd2_941b_3830_e000),
    (0x3ff5_4725_e2a7_7f93, 0xbfd2_3ec5_8cda_8800),
    (0x3ff5_2aff_4206_4583, 0xbfd1_e9e1_2927_9000),
    (0x3ff5_0f22_dbb2_bddf, 0xbfd1_956d_2b48_f800),
    (0x3ff4_f38f_4734_ded7, 0xbfd1_4167_9ab9_f800),
    (0x3ff4_d843_cfde_2840, 0xbfd0_edd0_94ef_9800),
    (0x3ff4_bd3e_c078_a3c8, 0xbfd0_9aa5_18db_1000),
    (0x3ff4_a27f_c3e0_258a, 0xbfd0_47e6_5263_b800),
    (0x3ff4_8805_24d4_8434, 0xbfcf_eb22_4586_f000),
    (0x3ff4_6dce_1b19_2d0b, 0xbfcf_474a_7517_b000),
    (0x3ff4_53d9_d339_1854, 0xbfce_a444_3d10_3000),
    (0x3ff4_3a27_44b4_845a, 0xbfce_020d_44e9_b000),
    (0x3ff4_20b5_4115_f8fb, 0xbfcd_60a2_2977_f000),
    (0x3ff4_0782_da3e_f4b1, 0xbfcc_c001_0495_9000),
    (0x3ff3_ee8f_5d57_fe8f, 0xbfcc_2029_5689_1000),
    (0x3ff3_d5d9_a00b_4ce9, 0xbfcb_8117_8d81_1000),
    (0x3ff3_bd60_c010_c12b, 0xbfca_e2c9_ccd3_d000),
    (0x3ff3_a524_2b75_dab8, 0xbfca_4540_2e12_9000),
    (0x3ff3_8d22_cd9f_d002, 0xbfc9_a877_681d_f000),
    (0x3ff3_755b_c584_7a1c, 0xbfc9_0c6d_6948_3000),
    (0x3ff3_5dce_49ad_36e2, 0xbfc8_7120_a645_c000),
    (0x3ff3_4679_984d_d440, 0xbfc7_d68f_b414_3000),
    (0x3ff3_2f5c_ceff_cb24, 0xbfc7_3cb8_3c62_7000),
    (0x3ff3_1877_75a1_0d49, 0xbfc6_a39a_9b37_6000),
    (0x3ff3_01c8_373e_3990, 0xbfc6_0b31_54b7_a000),
    (0x3ff2_eb4e_bb95_f841, 0xbfc5_737d_7624_3000),
    (0x3ff2_d50a_0219_a9d1, 0xbfc4_dc7b_8fc2_3000),
    (0x3ff2_bef9_a8b7_fd2a, 0xbfc4_462c_51d2_0000),
    (0x3ff2_a91c_7a0c_1bab, 0xbfc3_b08a_bc83_0000),
    (0x3ff2_9372_6014_b530, 0xbfc3_1b99_6b49_0000),
    (0x3ff2_7dfa_5757_a1f5, 0xbfc2_8754_90a4_4000),
    (0x3ff2_68b3_9b1d_3bbf, 0xbfc1_f3b9_f879_a000),
    (0x3ff2_539d_838f_f5bd, 0xbfc1_60c8_252c_a000),
    (0x3ff2_3eb7_aac9_083b, 0xbfc0_ce7f_57f7_2000),
    (0x3ff2_2a01_2ba9_40b6, 0xbfc0_3cdc_49fe_a000),
    (0x3ff2_1579_96cc_4132, 0xbfbf_57bd_bc4b_8000),
    (0x3ff2_0120_1dd2_fc9b, 0xbfbe_3708_9640_4000),
    (0x3ff1_ecf4_494d_480b, 0xbfbd_1798_3ef9_4000),
    (0x3ff1_d8f5_528f_6569, 0xbfbb_f967_4ed8_a000),
    (0x3ff1_c523_1157_7e7c, 0xbfba_dc79_202f_6000),
    (0x3ff1_b17c_74cb_26e9, 0xbfb9_c0c3_e728_8000),
    (0x3ff1_9e01_0c2c_1ab6, 0xbfb8_a646_b372_c000),
    (0x3ff1_8ab0_7bb6_70bd, 0xbfb7_8d01_b3ac_0000),
    (0x3ff1_778a_25ef_bcb6, 0xbfb6_74f1_4538_0000),
    (0x3ff1_648d_354c_31da, 0xbfb5_5e0e_6d87_8000),
    (0x3ff1_51b9_9027_5fdd, 0xbfb4_485c_dea1_e000),
    (0x3ff1_3f0e_a432_d24c, 0xbfb3_33d9_4d6a_a000),
    (0x3ff1_2c8b_7210_f9da, 0xbfb2_2079_f8c5_6000),
    (0x3ff1_1a30_28ec_b531, 0xbfb1_0e46_9862_2000),
    (0x3ff1_07fb_da84_34af, 0xbfaf_fa6c_6ad2_0000),
    (0x3ff0_f5ee_0f4e_6bb3, 0xbfad_da8d_4a77_4000),
    (0x3ff0_e406_5d2a_9fce, 0xbfab_bcec_e485_0000),
    (0x3ff0_d244_632c_a521, 0xbfa9_a189_4012_c000),
    (0x3ff0_c0a7_7ce2_981a, 0xbfa7_8858_3302_c000),
    (0x3ff0_af2f_83c6_36d1, 0xbfa5_715e_67d6_8000),
    (0x3ff0_9ddb_98a0_1339, 0xbfa3_5c8a_4965_8000),
    (0x3ff0_8cab_af52_e7df, 0xbfa1_49e3_6415_4000),
    (0x3ff0_7b9f_2f4e_28fb, 0xbf9e_72c0_82eb_8000),
    (0x3ff0_6ab5_8c35_8f19, 0xbf9a_55f1_5252_8000),
    (0x3ff0_59ee_a5ec_f92c, 0xbf96_3d62_cf81_8000),
    (0x3ff0_4949_cdd1_2c90, 0xbf92_28fb_8caa_0000),
    (0x3ff0_38c6_c6f0_ada9, 0xbf8c_317b_20f9_0000),
    (0x3ff0_2865_1379_32a9, 0xbf84_1935_5daa_0000),
    (0x3ff0_1824_27ea_7348, 0xbf78_1203_c2ec_0000),
    (0x3ff0_0804_0614_b195, 0xbf60_0409_7924_0000),
    (0x3fef_e01f_f726_fa1a, 0x3f6f_eff3_8490_0000),
    (0x3fef_a11c_c261_ea74, 0x3f87_dc41_353d_0000),
    (0x3fef_6310_b081_992e, 0x3f93_cea3_c4c2_8000),
    (0x3fef_25f6_3cee_adcd, 0x3f9b_9fc1_1489_0000),
    (0x3fee_e9c8_0391_13e7, 0x3fa1_b0d8_ce11_0000),
    (0x3fee_ae80_78cb_b1ab, 0x3fa5_8a5b_d001_c000),
    (0x3fee_741a_a29d_0c9b, 0x3fa9_5c83_40d8_8000),
    (0x3fee_3a91_830a_99b5, 0x3fad_276a_ef57_8000),
    (0x3fee_01e0_0960_9a56, 0x3fb0_7598_e598_c000),
    (0x3fed_ca01_e577_bb98, 0x3fb2_53f5_e30d_2000),
    (0x3fed_92f2_0b7c_9103, 0x3fb4_2edd_8b38_0000),
    (0x3fed_5cac_66fb_5cce, 0x3fb6_0659_8757_c000),
    (0x3fed_272c_aa5e_de9d, 0x3fb7_da76_356a_0000),
    (0x3fec_f26e_3e6b_2ccd, 0x3fb9_ab43_4e1c_6000),
    (0x3fec_be6d_a2a7_7902, 0x3fbb_78c7_bb0d_6000),
    (0x3fec_8b26_6d37_086d, 0x3fbd_4313_32e7_2000),
    (0x3fec_5894_bd5d_5804, 0x3fbf_0a31_71de_6000),
    (0x3fec_26b5_33bb_9f8c, 0x3fc0_6715_2b91_4000),
    (0x3feb_f583_eeec_e73f, 0x3fc1_4785_8292_b000),
    (0x3feb_c4fd_75db_96c1, 0x3fc2_266e_cdca_3000),
    (0x3feb_951e_0c86_4a28, 0x3fc3_03d7_a6c5_5000),
    (0x3feb_65e2_c5ef_3e2c, 0x3fc3_dfc3_3c33_1000),
    (0x3feb_3748_67c9_888b, 0x3fc4_ba36_6b7a_8000),
    (0x3feb_094b_211d_304a, 0x3fc5_9339_28d1_f000),
    (0x3fea_dbe8_85f2_ef7e, 0x3fc6_6acd_2418_f000),
    (0x3fea_af1d_3160_3da2, 0x3fc7_40f8_ec66_9000),
    (0x3fea_82e6_3fd3_58a7, 0x3fc8_15c0_f51a_f000),
    (0x3fea_5740_ef09_738b, 0x3fc8_e929_54f6_8000),
    (0x3fea_2c2a_90ab_4b27, 0x3fc9_bb36_02f8_4000),
    (0x3fea_01a0_1393_f2d1, 0x3fca_8bed_1c2c_0000),
    (0x3fe9_d79f_24db_3c1b, 0x3fcb_5b51_5c01_d000),
    (0x3fe9_ae25_05c7_b190, 0x3fcc_2967_ccbc_c000),
    (0x3fe9_852e_f297_ce2f, 0x3fcc_f635_d548_6000),
    (0x3fe9_5cba_eea4_4b75, 0x3fcd_c1bd_3446_c000),
    (0x3fe9_34c6_9de7_4838, 0x3fce_8c01_b8cf_e000),
    (0x3fe9_0d4f_2f67_52e6, 0x3fcf_5509_c017_9000),
    (0x3fe8_e652_8eff_d79d, 0x3fd0_0e6c_121f_b800),
    (0x3fe8_bfce_9fcc_007c, 0x3fd0_71b8_0e93_d000),
    (0x3fe8_99c0_dabe_c30e, 0x3fd0_d46b_9e86_7000),
    (0x3fe8_7427_aa23_17fb, 0x3fd1_3687_334b_d000),
    (0x3fe8_4f00_acb3_9a08, 0x3fd1_980d_6723_4800),
    (0x3fe8_2a49_e865_3e55, 0x3fd1_f8ff_e0cc_8000),
    (0x3fe8_0601_95f4_0260, 0x3fd2_595f_d763_6800),
    (0x3fe7_e225_63e0_a329, 0x3fd2_b930_0914_a800),
    (0x3fe7_beb3_77dc_b5ad, 0x3fd3_1872_1043_6000),
    (0x3fe7_9baa_6797_25c2, 0x3fd3_7726_6dec_1800),
    (0x3fe7_7907_f217_0657, 0x3fd3_d54f_fbaf_3000),
    (0x3fe7_56ca_dbd6_130c, 0x3fd4_32ee_e32f_e000),
];

/// `__log_data.poly`: `A[0..5]`.
const LOG_POLY: [u64; 5] = [
    0xbfe0_0000_0000_0001,
    0x3fd5_5555_5551_305b,
    0xbfcf_ffff_ffeb_4590,
    0x3fc9_99b3_24f1_0111,
    0xbfc5_5575_e506_c89f,
];

/// `__log_data.poly1`: `B[0..11]`.
const LOG_POLY1: [u64; 11] = [
    0xbfe0_0000_0000_0000,
    0x3fd5_5555_5555_5577,
    0xbfcf_ffff_ffff_fdcb,
    0x3fc9_9999_9995_dd0c,
    0xbfc5_5555_5567_45a7,
    0x3fc2_4924_a344_de30,
    0xbfbf_ffff_a442_3d65,
    0x3fbc_7184_282a_d6ca,
    0xbfb9_99eb_43b0_68ff,
    0x3fb7_8182_f7af_d085,
    0xbfb5_5213_75d1_45cd,
];

/// `__exp_data.poly`: `C2..C5`.
const EXP_POLY: [u64; 4] = [
    0x3fdf_ffff_ffff_fdbd,
    0x3fc5_5555_5555_543c,
    0x3fa5_5555_cf17_2b91,
    0x3f81_1111_67a4_d017,
];
const EXP_INVLN2N: u64 = 0x4067_1547_652b_82fe;
const EXP_SHIFT: u64 = 0x4338_0000_0000_0000;
const EXP_NEGLN2HIN: u64 = 0xbf76_2e42_fefa_0000;
const EXP_NEGLN2LON: u64 = 0xbd0c_f79a_bc9e_3b3a;
const LOG_LN2HI: u64 = 0x3fe6_2e42_fefa_3800;
const LOG_LN2LO: u64 = 0x3d2e_f357_93c7_6730;

/// `1.0`.
const ONE: f64 = 1.0;

/// `LO` of the near-1 test of `log`: the bits of `1.0 - 0x1p-4`.
const LOG_NEAR_ONE_LO: u64 = 0x3fee_0000_0000_0000;

/// `HI - LO` of the near-1 test of `log`: `1.0 + 0x1.09p-4` minus `LO`.
const LOG_NEAR_ONE_SPAN: u64 = 0x0003_0900_0000_0000;

/// `OFF` of `log`: the bits of `0x1.6p-1`.
const LOG_OFF: u64 = 0x3fe6_0000_0000_0000;

/// `0x1p27`, the splitting factor of `log` near 1.
const TWO_POW_27: f64 = 134_217_728.0;

/// `0x1p52`, the normalisation factor of a subnormal `log` argument.
const TWO_POW_52: f64 = 4_503_599_627_370_496.0;

/// `0x1p1009`: the scale of an `exp` result that may overflow.
const TWO_POW_1009: u64 = 0x7f00_0000_0000_0000;

/// `0x1p-1022`: the scale of an `exp` result that may underflow.
const TWO_POW_MINUS_1022: u64 = 0x0010_0000_0000_0000;

/// `__math_oflow(0)`: `0x1p769 * 0x1p769`, positive infinity.
fn overflow() -> f64 {
    let y = f64::from_bits(0x7000_0000_0000_0000);
    y * y
}

/// `__math_uflow(0)`: `0x1p-767 * 0x1p-767`, positive zero.
fn underflow() -> f64 {
    let y = f64::from_bits(0x1000_0000_0000_0000);
    y * y
}

/// `exp`'s `specialcase`: `scale * (1 + tmp)` for a reduction index `ki`
/// whose scale lies outside the normal range, as `__ieee754_exp_fma` computes
/// it (`+256` to `+476` of the function).
fn exp_special_case(tmp: f64, sbits: u64, ki: u64) -> f64 {
    if ki & 0x8000_0000 == 0 {
        // k > 0: `0x1p1009 * fma(scale, tmp, scale)`.
        let scale = f64::from_bits(sbits.wrapping_sub(1009 << 52));
        return scale.mul_add(tmp, scale) * f64::from_bits(TWO_POW_1009);
    }
    // k < 0: `scale * tmp` is used twice and not fused.
    let scale = f64::from_bits(sbits.wrapping_add(1022 << 52));
    let product = tmp * scale;
    let y = scale + product;
    // `comisd`/`ja`: the branch is taken only for `1 > y` (never for a NaN).
    if ONE > y {
        let hi = y + ONE;
        let lo = (scale - y) + product;
        let lo = ((ONE - hi) + y) + lo;
        let y = (lo + hi) - ONE;
        if y == 0.0 {
            return 0.0;
        }
        return y * f64::from_bits(TWO_POW_MINUS_1022);
    }
    y * f64::from_bits(TWO_POW_MINUS_1022)
}

/// `e` to the power `x` as `__ieee754_exp_fma` of the reference build's GNU C
/// Library 2.39 computes it, NaN bits included (module documentation).
pub(crate) fn exp(x: f64) -> f64 {
    let ix = x.to_bits();
    let mut abstop = ((ix >> 52) & 0x7ff) as u32;
    // `abstop - top12(0x1p-54) >= top12(512) - top12(0x1p-54)`, unsigned.
    if abstop.wrapping_sub(0x3c9) > 0x3e {
        if abstop < 0x3c9 {
            // |x| < 2^-54 (0 included): `x + 1.0`.
            return x86_64::add(x, ONE);
        }
        if abstop >= 0x409 {
            if ix == 0xfff0_0000_0000_0000 {
                return 0.0;
            }
            if abstop == 0x7ff {
                return x86_64::add(x, ONE);
            }
            return if ix >> 63 != 0 {
                underflow()
            } else {
                overflow()
            };
        }
        // 512 <= |x| < 1024: the special case below.
        abstop = 0;
    }
    let shift = f64::from_bits(EXP_SHIFT);
    // `kd = fma(x, N/ln2, shift)`: the rounding shift is fused.
    let kd = x.mul_add(f64::from_bits(EXP_INVLN2N), shift);
    let ki = kd.to_bits();
    let kd = kd - shift;
    // `r = fma(kd, -ln2lo/N, fma(kd, -ln2hi/N, x))`.
    let r = kd.mul_add(
        f64::from_bits(EXP_NEGLN2LON),
        kd.mul_add(f64::from_bits(EXP_NEGLN2HIN), x),
    );
    let index = 2 * (ki & 0x7f) as usize;
    let top = ki << 45;
    let tail = f64::from_bits(EXP_TAB[index]);
    let sbits = EXP_TAB[index + 1].wrapping_add(top);
    let c = EXP_POLY.map(f64::from_bits);
    // `tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5)`, fused
    // as `fma(r2 * r2, fma(r, C5, C4), fma(fma(r, C3, C2), r2, r + tail))`.
    let low = r.mul_add(c[1], c[0]);
    let head = r + tail;
    let r2 = r * r;
    let high = r.mul_add(c[3], c[2]);
    let sum = low.mul_add(r2, head);
    let r4 = r2 * r2;
    let tmp = r4.mul_add(high, sum);
    if abstop == 0 {
        return exp_special_case(tmp, sbits, ki);
    }
    let scale = f64::from_bits(sbits);
    scale.mul_add(tmp, scale)
}

/// `__math_invalid`: `(x - x) / (x - x)` with SSE's NaN rule.
fn invalid(x: f64) -> f64 {
    let difference = x86_64::sub(x, x);
    x86_64::div(difference, difference)
}

/// The natural logarithm of `x` as `__ieee754_log_fma` of the reference
/// build's GNU C Library 2.39 computes it, NaN bits included (module
/// documentation).
pub(crate) fn log(x: f64) -> f64 {
    let mut ix = x.to_bits();
    if ix.wrapping_sub(LOG_NEAR_ONE_LO) < LOG_NEAR_ONE_SPAN {
        return log_near_one(x, ix);
    }
    let top = (ix >> 48) as u32;
    if top.wrapping_sub(0x0010) > 0x7fdf {
        // x < 0x1p-1022, infinite or NaN.
        if ix.wrapping_mul(2) == 0 {
            // `__math_divzero(1)`: `-1.0 / 0.0`.
            return f64::NEG_INFINITY;
        }
        if ix == 0x7ff0_0000_0000_0000 {
            return x;
        }
        if top & 0x8000 != 0 || top & 0x7ff0 == 0x7ff0 {
            return invalid(x);
        }
        // A positive subnormal, normalised.
        ix = (x * TWO_POW_52).to_bits().wrapping_sub(52 << 52);
    }
    let tmp = ix.wrapping_sub(LOG_OFF);
    let i = ((tmp >> 45) & 0x7f) as usize;
    let k = (tmp as i64) >> 52;
    let iz = ix.wrapping_sub(tmp & (0xfff << 52));
    let (invc, logc) = LOG_TAB[i];
    let z = f64::from_bits(iz);
    let a = LOG_POLY.map(f64::from_bits);
    // `vcvtsi2sd` of the low 32 bits of `k`.
    let kd = f64::from(k as i32);
    // `w = fma(kd, ln2hi, logc)`, `r = fma(z, invc, -1)`.
    let w = kd.mul_add(f64::from_bits(LOG_LN2HI), f64::from_bits(logc));
    let r = z.mul_add(f64::from_bits(invc), -ONE);
    let q1 = r.mul_add(a[2], a[1]);
    let hi = r + w;
    let r2 = r * r;
    // `lo = fma(kd, ln2lo, w - hi + r)`.
    let lo = kd.mul_add(f64::from_bits(LOG_LN2LO), (w - hi) + r);
    let r3 = r * r2;
    let q2 = r.mul_add(a[4], a[3]);
    let lo = r2.mul_add(a[0], lo);
    let q = q2.mul_add(r2, q1);
    // `y = lo + r2 * A[0] + r * r2 * (A[1] + r * A[2] + r2 * (A[3] + r * A[4])) + hi`.
    r3.mul_add(q, lo) + hi
}

/// `log` for `x` in `[1 - 0x1p-4, 1 + 0x1.09p-4)`, as `__ieee754_log_fma`
/// computes it (`+256` to `+466` of the function).
fn log_near_one(x: f64, ix: u64) -> f64 {
    if ix == ONE.to_bits() {
        return 0.0;
    }
    let b = LOG_POLY1.map(f64::from_bits);
    let r = x - ONE;
    let p1 = r.mul_add(b[2], b[1]);
    let p4 = r.mul_add(b[5], b[4]);
    let r2 = r * r;
    let p7 = r.mul_add(b[8], b[7]);
    let p1 = r2.mul_add(b[3], p1);
    let p4 = r2.mul_add(b[6], p4);
    let r3 = r * r2;
    let p7 = r2.mul_add(b[9], p7);
    let p7 = r3.mul_add(b[10], p7);
    let p4 = p7.mul_add(r3, p4);
    // `y = r3 * (B[1] + r * B[2] + r2 * B[3] + r3 * (B[4] + ...))`, added to
    // `lo` below in one fused step.
    let poly = p4.mul_add(r3, p1);
    // `w = r * 0x1p27; rhi = r + w - w`, both uses of `w` fused:
    // `rhi = fnma(0x1p27, r, fma(r, 0x1p27, r))`.
    let sum = r.mul_add(TWO_POW_27, r);
    let rhi = (-TWO_POW_27).mul_add(r, sum);
    let square = rhi * rhi;
    let rlo = r - rhi;
    // `w = rhi * rhi * B[0]; hi = r + w; lo = r - hi + w`, each use of `w`
    // fused.
    let hi = square.mul_add(b[0], r);
    let lo = square.mul_add(b[0], r - hi);
    // `lo += B[0] * rlo * (rhi + r)`.
    let lo = (b[0] * rlo).mul_add(r + rhi, lo);
    poly.mul_add(r3, lo) + hi
}

/// The arc tangent of `x`: the host's on x86_64 Linux with the GNU C Library
/// (module documentation: the reference build's `__atan_fma` has no
/// licence-clean upstream; exact where the host selects it).
#[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
pub(crate) fn atan(x: f64) -> f64 {
    x.atan()
}

/// The arc tangent of `x`: the `libm` crate's where the host is not x86_64
/// Linux with the GNU C Library, with `__atan_fma`'s `x + x` for a NaN (module
/// documentation).
#[cfg(not(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64")))]
pub(crate) fn atan(x: f64) -> f64 {
    if x.is_nan() {
        x86_64::add(x, x)
    } else {
        libm::atan(x)
    }
}

/// `sqrt`, with the NaN bits of `sqrtsd`.
pub(crate) fn sqrt(x: f64) -> f64 {
    x86_64::sqrt(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> String {
        std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/feature_finder_picked/glibc_libm_probe.tsv"
        ))
        .unwrap()
    }

    fn hex(text: &str) -> u64 {
        u64::from_str_radix(text, 16).unwrap()
    }

    /// Whether [`atan`] is the host's and the host can select `__atan_fma`
    /// (x86_64 Linux with the GNU C Library, on a CPU with FMA; module
    /// documentation). Elsewhere [`atan`] is the documented fallback, which
    /// is not the reference.
    ///
    /// This asks about the **processor**, because glibc's indirect-function
    /// resolver reads CPUID when `libm.so.6` is loaded. It must therefore not
    /// be asked with `std::arch::is_x86_feature_detected!("fma")`, which is
    /// documented to answer `true` *without* consulting the processor whenever
    /// the feature is already enabled at compile time — and since
    /// `.cargo/config.toml` builds x86_64 with `-C target-feature=+fma`, that
    /// is now every build of this crate, so the macro would fold to a constant
    /// here (`docs/FMA_BUILD_FLAG.md` section 6). It asks `raw_cpuid` instead,
    /// which always executes `cpuid`, so the question stays the one that was
    /// asked.
    ///
    /// **Why not `system::cpu_features::cpu_provides_fma`**, which is the
    /// production copy of exactly this read: `analysis` may not name `crate::
    /// system`. `system` already reaches `analysis` through `format`, so that
    /// edge closes a module cycle and `tools/check_module_cycles.py` fails on
    /// it. The two must stay the same architectural bit — leaf 1, `ECX` bit 12,
    /// absent leaf 1 counting as present — and `cpu_features`' own
    /// documentation is where that convention is written down.
    ///
    /// Nothing about what the test compares changes: a `+fma` binary cannot
    /// start on a processor without FMA, so the two mechanisms agree on every
    /// host that can run this code at all. Only the mechanism differs, and this
    /// one cannot be silenced by a build flag.
    fn atan_is_reference() -> bool {
        #[cfg(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64"))]
        {
            raw_cpuid::CpuId::new()
                .get_feature_info()
                .is_none_or(|info| info.has_fma())
        }
        #[cfg(not(all(target_os = "linux", target_env = "gnu", target_arch = "x86_64")))]
        {
            false
        }
    }

    fn splitmix64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn fnv64(mut h: u64, v: u64) -> u64 {
        for i in 0..8 {
            h ^= (v >> (8 * i)) & 0xff;
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        h
    }

    /// `(d >> 11) * 2^-53`.
    fn unit(d: u64) -> f64 {
        (d >> 11) as f64 * f64::from_bits(0x3ca0_0000_0000_0000)
    }

    fn scaled(d: u64, span: u64, offset: i64) -> f64 {
        let mut mix = d;
        let e = (((d >> 52) % span) as i64 + offset + 1023) as u64;
        let sign = splitmix64(&mut mix) >> 63;
        f64::from_bits((sign << 63) | (e << 52) | (d & 0x000f_ffff_ffff_ffff))
    }

    /// The generators of `libm_probe.c` (`../oracle/ffap-complete-fix3`).
    fn draw(set: u32, state: &mut u64) -> f64 {
        let d = splitmix64(state);
        match set {
            11 | 23 | 31 => f64::from_bits(d),
            12 => -760.0 + unit(d) * 1480.0,
            13 => scaled(d, 71, -60),
            14 => -unit(d) * 50.0,
            15 => (if d & 1 != 0 { -1.0 } else { 1.0 }) * (512.0 + unit(d) * 512.0),
            21 => f64::from_bits(d & 0x7fff_ffff_ffff_ffff),
            22 => 1.0 + (unit(d) - 0.5) * 0.25,
            24 => ((d >> 11) + 1) as f64 * f64::from_bits(0x3ca0_0000_0000_0000),
            25 => f64::from_bits(d & 0x000f_ffff_ffff_ffff),
            32 => unit(d) * 100.0,
            _ => scaled(d, 81, -40),
        }
    }

    fn function(set: u32) -> fn(f64) -> f64 {
        if set < 20 {
            exp
        } else if set < 30 {
            log
        } else {
            atan
        }
    }

    #[test]
    fn the_fixture_is_the_reference_library() {
        let text = fixture();
        assert!(text.contains("where\tglibc\t2.39\n"));
        assert!(text.contains("where\texp\t0x3a700\texp\t/lib/x86_64-linux-gnu/libm.so.6\n"));
        assert!(text.contains("where\tlog\t0x3a5f0\tlogf64\t"));
        assert!(text.contains("where\tatan\t0x7a860\t(none)\t"));
    }

    /// Every special value through `exp`, `log` and `atan`: bit for bit, NaN
    /// bits included (`atan` only where it is the reference's, see
    /// `atan_is_reference`).
    #[test]
    fn special_values_match_the_executed_library() {
        let mut compared = 0;
        for line in fixture().lines().filter(|line| line.starts_with("grid\t")) {
            let fields: Vec<&str> = line.split('\t').collect();
            let x = f64::from_bits(hex(fields[2]));
            let expected = hex(fields[3]);
            let actual = match fields[1] {
                "exp" => exp(x),
                "log" => log(x),
                "atan" if atan_is_reference() || x.is_nan() => atan(x),
                "atan" => continue,
                other => panic!("{other}"),
            };
            assert_eq!(
                actual.to_bits(),
                expected,
                "{}({:016x}): {actual:e}",
                fields[1],
                x.to_bits()
            );
            compared += 1;
        }
        assert!(compared >= 160, "{compared}");
    }

    /// The first 256 inputs of every generated set, and the digest of the
    /// first `2^20`.
    #[test]
    fn generated_inputs_match_the_executed_library() {
        let text = fixture();
        let mut rows: std::collections::BTreeMap<u32, Vec<(u64, u64)>> =
            std::collections::BTreeMap::new();
        let mut digests = std::collections::BTreeMap::new();
        for line in text.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields[0] {
                "row" => rows
                    .entry(fields[1].parse().unwrap())
                    .or_default()
                    .push((hex(fields[2]), hex(fields[3]))),
                "digest" if fields[2] == "0" => {
                    digests.insert(fields[1].parse::<u32>().unwrap(), hex(fields[3]));
                }
                _ => {}
            }
        }
        assert_eq!(rows.len(), 13);
        for (&set, expected) in &rows {
            if set >= 30 && !atan_is_reference() {
                continue;
            }
            let f = function(set);
            let mut state = u64::from(set);
            let mut digest = 0xcbf2_9ce4_8422_2325u64;
            for i in 0..(1usize << 20) {
                let x = draw(set, &mut state);
                let y = f(x);
                if let Some(&(ex, ey)) = expected.get(i) {
                    assert_eq!(x.to_bits(), ex, "set {set} input {i}");
                    assert_eq!(y.to_bits(), ey, "set {set} input {i}: {x:e}");
                }
                digest = fnv64(fnv64(digest, x.to_bits()), y.to_bits());
            }
            assert_eq!(digest, digests[&set], "set {set}");
        }
    }

    /// The disassembled special cases, by hand: `exp` of `-inf`, of huge and
    /// tiny arguments, and `log` of zero, one, infinity and negative numbers.
    #[test]
    fn special_cases_follow_the_disassembly() {
        assert_eq!(exp(f64::NEG_INFINITY).to_bits(), 0);
        assert_eq!(exp(f64::INFINITY), f64::INFINITY);
        assert_eq!(exp(1000.0), f64::INFINITY);
        assert_eq!(exp(-1000.0).to_bits(), 0);
        assert_eq!(exp(0.0), 1.0);
        assert_eq!(exp(-0.0), 1.0);
        assert_eq!(log(1.0).to_bits(), 0);
        assert_eq!(log(0.0), f64::NEG_INFINITY);
        assert_eq!(log(-0.0), f64::NEG_INFINITY);
        assert_eq!(log(f64::INFINITY), f64::INFINITY);
        // `(x - x) / (x - x)`: SSE's default NaN for a finite negative `x` and
        // for `-inf`, the quieted `x` for a NaN.
        assert_eq!(log(-1.0).to_bits(), 0xfff8_0000_0000_0000);
        assert_eq!(log(f64::NEG_INFINITY).to_bits(), 0xfff8_0000_0000_0000);
        assert_eq!(
            log(f64::from_bits(0x7ff0_0000_0000_0001)).to_bits(),
            0x7ff8_0000_0000_0001
        );
        assert_eq!(sqrt(-1.0).to_bits(), 0xfff8_0000_0000_0000);
        assert_eq!(
            atan(f64::from_bits(0xfff4_0000_0000_0000)).to_bits(),
            0xfffc_0000_0000_0000
        );
    }
}
