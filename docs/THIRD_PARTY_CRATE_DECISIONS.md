# Third-party dependency decisions

Where the C++ takes code from a third-party library, the port looks for a
suitable crate before writing its own (see `docs/PORTING_STATUS.md`). This
register records every such decision so it is not re-litigated by the next
porting agent, and states what would reopen it.

A crate is suitable only if it passes all five tests:

1. **MSRV**: builds under Rust 1.85, measured, not read off `rust-version`.
2. **Pure Rust**: no C/C++/Fortran toolchain in the build.
3. **Maintained**: recent releases and real downstream use.
4. **Fidelity**: reproduces the upstream class-test expectations at their
   asserted tolerances. Where a C++ oracle was compiled, the crate must not be
   further from it than the code it replaces.
5. **Determinism**: no runtime CPU dispatch on the executed path; same bits on
   aarch64 macOS and x86_64 Linux except where the platform `libm` already makes
   the existing code differ.

Evidence comes from a survey run on 2026-09-13: three independent finders and
one evaluation per candidate. Every evaluation measured MSRV on the Linux build
node and fidelity in a throwaway crate, often against an executed Boost or Eigen
oracle. The survey changed nothing in the repository. "Pending" means the
decision is made but the refactor has not landed. Crate wave 1 (2026-09-14)
landed three refactors, part of a fourth, and overturned the digamma decision;
each row names the commit that carries it.

## Replace with a crate

| C++ library and use | Rust today | Crate, pin, features | Status | Measured basis |
|---|---|---|---|---|
| evergreen FFT (`KernelDensityEstimation`) | `src/math/fft.rs` | `rustfft =6.4.1`, always `FftPlannerScalar` | Done, `d396f42` | Upstream KDE expectations pass unchanged. The SIMD planner differs from scalar in 92-99% of components for at most 1.1x on KDE lengths. |
| Eigen `unsupported/NonLinearOptimization` Levenberg-Marquardt (Gauss, Gamma, Gumbel, Gumbel-MLE fitters; later `TraceFitter`, `LevMarqFitter1D`) | `src/math/fitters/levenberg_marquardt.rs` | `levenberg-marquardt =0.14.0`; `nalgebra =0.33.3`, `default-features = false`, `alloc` + `libm` | Pending | The swapped solver passes all 30 fitter, 9 PEP and 28 unit tests under 1.85. On the one C++-published Gauss case it is slightly closer to C++ than the port (A 1.5e-12 vs 2.2e-12). Bit-identical across machines with `libm` residuals; `matrixmultiply` stays out of the graph. |
| Boost.Math normal quantile (`MultipleTesting` lfdr) | `standard_normal_quantile` in `src/math/multiple_testing.rs` | `statrs =0.18.0`, `default-features = false`; call `erfc_inv` directly | Done, `823c8e9` (merge of `29bfbc4`) | `-SQRT_2 * erfc_inv(2p)` is Boost's own formula, now in Boost's statement order. Measured after the refactor: at most 1.55 ulp from correctly rounded quantiles over the 40 reference points (2.93e-16 relative), where the replaced Acklam+Halley was up to 1.1e-9 relative off; bit-identical to Boost 1.92's double quantile at 79.7% of 97,489 grid points on macOS arm64, against 46.1% before. The survey had measured 4.4e-16 against 1.1e-9 relative to Boost 1.92. lfdr tests unchanged. |
| Boost.Random `mt19937_64` engine (`DecoyGenerator` via `Math::RandomShuffler`; `std::mt19937_64` in `UniqueIdGenerator`) | engine in `src/chemistry/decoy_random.rs` | `rand_mt =6.0.3`, `default-features = false` | Done, `f913cf5` (merge of `1d0e9c8`) | 0 differing words in 52M across 20,012 seeds. The 13 upstream DecoyGenerator strings pass. Boost's `normalize` never affects output. Keep the Boost `uniform_int` bucket mapping and the Fisher-Yates loop: `rand`'s mapping changes every decoy. Refactor check: 0 of 33,310,000 words and 0 of 249,750 shuffles differ from the replaced engine; the engine size stays 2504 bytes. |
| Boost.Regex (`SpectrumLookup`, `SpectrumNativeIDParser`, `MascotXMLFile`, `PepXMLFile`, `MzIdentMLHandler`, `MzTabFile`, `IndexedMzMLDecoder`, `FalseDiscoveryRate`, `EnzymaticDigestion`, `RNaseDigestion`) | Hand-written per-pattern scanners; arbitrary expressions refused with `Unsupported` | `fancy-regex =0.19.2`, behind one Boost-compatible facade | Pending; the facade is in review on `crate/regex-facade` and not merged | The facade needs ASCII bytes mode, dot-matches-newline, `^`/`$` rewritten to Boost's line anchors (`\n`, `\r`, `\f`) and a `match_not_initial_null` retry. With those, 0 mismatches on 240,223 cases over 65 source patterns. The same corpus exposes three port bugs: `mzidentml` `split_annotation` drops multi-loss annotations, `mztab` misses nested brackets, `pepxml` anchors `$` at end of string and parses u64 where C++ uses int32. The plain `regex` crate rejects the lookaround enzyme expressions and duplicate group names. |
| Xerces character/entity reference expansion on read | `expand_reference` in `mascot_xml.rs`; inline copies in `identification_xml.rs` and `mzidentml.rs` | `quick-xml` (existing): `BytesRef::resolve_char_ref`, `escape::resolve_xml_entity` | Partly done, `db19767` (merge of `be2bb15`): `identification_xml.rs` only | 15/20 probes identical. The 5 that differ (`&#X2E;`, `&#+46;`, `&#x+2E;`, `&#0;`, `&#x0;`) are accepted by the hand-rolled code and rejected by quick-xml and libxml2, as the XML 1.0 grammar requires. In `identification_xml.rs` element text now resolves through `resolve_char_ref` and `resolve_xml_entity`; the only accepted-set change there was signed references (`&#+46;`, `&#x+2E;`), because `&#X2E;`, `&#0;` and `&#x0;` were already refused (14,141,319 reference names compared). Attribute values still go through `escape::unescape`, which uses the same character-reference parser. Remaining: `expand_reference` in `mascot_xml.rs` and the inline copy in `mzidentml.rs`. |
| libcurl URL parsing (`NetworkGetRequest`) | host-stripping in `check_http_url` | `ureq::http::Uri` (existing, via ureq) | Done, `db19767` (merge of `be2bb15`) | 9/9 existing cases pass. Over 2M random URLs the new check never rejects a usable host, and it closes 275 empty-host URLs that reached DNS. Keep the scheme gate: `http::Uri` rejects `file://` as malformed. Measured after the refactor on the same 2M URLs: 0 previously refused URLs are accepted, the 275 empty-host URLs fail offline, and 866,415 unparseable URLs are refused by the preflight instead of by ureq, with the same error variant. |

## Keep the port's own code

| C++ library and use | Rust | Why no crate |
|---|---|---|
| eol-bspline (Ooyama smoothing B-spline, `BSpline2d`) | `src/processing/spline/b_spline.rs` | No crate implements Ooyama's penalised uniform-node spline. `csaps` and `splinefit` fail all 41 `TransformationModelBSpline` evaluations; the port matches the C++ probe bit for bit. |
| IsoSpec (fine isotope patterns) | `src/chemistry/fine_isotopes.rs` | No IsoSpec port or binding exists. `mzcore`/`rustyms`/`chemical_elements`/`mscore` compute coarse or binomial distributions; `mscore` fails all 14 fructose rows at 1e-7. |
| Eigen `JacobiSVD` (Savitzky-Golay coefficients) | `src/processing/smoothing.rs` | `nalgebra` and `faer` are further from Eigen (to 2.2e-10 and 8.6e-11 vs 1.6e-12) and give different bits on different machines. |
| LIBSVM prediction (FeatureFindingMetabo isotope models) | `src/analysis/feature_finding_metabo/predictor.rs` | `ffsvm` computes in f32: 9-11 of 488 executed-LIBSVM labels flip, differently per machine. The other crates cannot load LIBSVM models or are abandoned C bindings. |
| Eigen `SparseVector<float>` (`BinnedSpectrum`) | `src/comparison.rs` | `nalgebra-sparse` has no sparse vector. `sprs` cannot accumulate at random indices, so a map would stay; 0.11.5 needs Rust 1.88. |
| GTE `ApprHeightLine2` (linear transformation model) | `LinearModel::fit` in `src/analysis/transformations.rs` | `linreg` matches bit for bit but last released in 2019 and pulls in syn 1. `linregress` uses an SVD and misses an exact upstream assertion. |
| Boost.Math binomial complement (`MathFunctions`) | `binomial_cdf_complement` in `src/concept/math_functions.rs` | `statrs` caps its continued fraction at 140 steps and silently returns wrong values (-0.54 at N=1e8). `special` is no closer to Boost. |
| Boost.Math digamma (`GammaDistributionFitter` Jacobian) | `digamma` in `src/math/fitters/gamma.rs` | `special =0.14.1` (`no_std`) was tried on `crate/digamma` (`d94274c`) and dropped: it fails the fidelity test against the executed libOpenMS. The survey's 4/3 against 11/7 ulp were measured against a rebuilt Eigen 3.5.0 + Boost 1.92 probe, which is itself +7/+5 ulp from the shipped libOpenMS. Against the executed libOpenMS from the class-test start (1, 3), the crate lands -8/-5 ulp (b/p) on macOS arm64 and +4/+3 on Linux x86_64, while the hand-rolled series lands -4/-2 on macOS. On the oracle's own platform the crate is therefore further on both parameters, and on Linux it is 1 ulp further in p. From the default start (1, 5) the crate is +1/+1 on both machines, closer than the series. With the crate, macOS and Linux also differ by 12/8 ulp from (1, 3). Enabling `special/std` together with `no_std` recurses until the stack overflows. Reopen if a crate is no further from the executed libOpenMS than the series from every class-test start. |
| Boost.Math normal pdf (`GaussFitter`, `SpectrumCheapDPCorr`) | `src/math/fitters/gauss.rs`, `src/comparison.rs` | The `gauss.rs` transcription is bit-identical to Boost on 2M inputs; `statrs` is not. `comparison.rs` used a `ROOT_TWO_PI` literal one ulp above Boost's run-time `sqrt(2*pi)`, a port bug rather than a crate question. Fixed in `29bfbc4` (now `SQRT_TWO_PI = (2.0 * PI).sqrt()`), which moves the SpectrumCheapDPCorr cross score by 2 ulp. |
| Xerces encoding detection and XML well-formedness | readers across `src/format` | quick-xml's `encoding` feature removes `unescape_value` and decodes ISO-8859-1 as windows-1252. `xmlparser` and `roxmltree` miss Xerces rejections. |
| `XMLHandler::writeXMLEscape` and friends (write side) | escapers in `paramxml`, `identification_xml`, `qcml`, `pepxml`, ... | OpenMS's own code, not a library. Its byte formats (`&#x9;`, no `&apos;`) are asserted by tests, and escaping is metered per byte. |
| libc `fnmatch` / `PathMatchSpecA` (`File::fileList`) | `src/system/file.rs` | A system library, and platform-inconsistent (macOS and glibc disagree on 117 of 300k pairs). `glob` and `globset` fail the port's tests (no escapes, byte matching, `{}` alternation). |
| libc `timegm`/`gmtime_r` calendar arithmetic | `src/data_structures/datetime.rs` | `chrono`'s `NaiveDate` stops at year 262143; the C++ probe covers the full i32 range. Only the day-validation tables in `acquisition.rs`/`log_stream.rs` and the date in `get_unique_name` can move to `chrono`, with no new dependency. The `acquisition.rs` and `log_stream.rs` day checks now use `NaiveDate::from_ymd_opt` (`db19767`, merge of `be2bb15`), with 0 disagreements against the replaced tables over 655,360,000 year, month and day combinations; `get_unique_name` in `file.rs` is untouched. |

## Combined dependency check

The six pending crates were added together to a copy of the tree and checked on
the Linux build node on 2026-09-13, before any refactor. The pins are exactly
those in the table above. Results:

- `cargo +1.85.0 check` succeeds with `--all-features --all-targets` and with
  `--no-default-features`.
- Feature unification resolves `nalgebra` to `alloc` and `libm` only, and
  `special` to `no_std` only.
- `statrs` brings in no `nalgebra`.
- None of `matrixmultiply`, `wide` or `safe_arch`, the crates that pick CPU
  kernels at runtime, is in the graph.
- The committed `Cargo.lock` gains 19 packages. Five of them (`ellip`,
  `lambert_w`, `num-lazy`, `numeric_literals` and `syn` 1) are optional
  dependencies behind features this crate does not enable. The lockfile records
  them, but they are never compiled: `cargo tree --target all -e all -i` finds
  no path to any of them.
- With that lockfile, `cargo +1.85.0 check --locked` passes with all features
  and with no default features, and the feature-boundary gate passes.

`special` has since been removed. The digamma decision keeps the port's own
series, so the crate wave 1 integration took `special` out of `Cargo.toml`, and
`Cargo.lock` lost `special` and the five packages only it required (`ellip`,
`lambert_w`, `num-lazy`, `numeric_literals` and `syn` 1). The other five pins
are unchanged.

## Constraints that would reopen a decision

- **MSRV.** Rust 1.85 forces older pins. A decision to raise the MSRV would
  unlock these updates, each needing its measurement re-run:
  - levenberg-marquardt 0.15 needs nalgebra 0.34 (Rust 1.87);
  - statrs 0.19 needs Rust 1.89;
  - quick-xml 0.42 needs Rust 1.86.
- **Evaluation budget.** The levenberg-marquardt crate counts only trust-region
  evaluations (`patience * (n + 1)`). Eigen's `maxfev` also counts the `n + 1`
  evaluations of `NumericalDiff` per Jacobian. Only fits that exhaust their
  budget can differ. The adapter must reproduce Eigen's accounting before the
  FEATUREFINDER trace fitters rely on it.
- **Feature unification.** `nalgebra/std` enabled anywhere in the graph silently
  brings back `matrixmultiply`'s runtime dispatch. Check `cargo tree -e features`
  when adding dependencies.
