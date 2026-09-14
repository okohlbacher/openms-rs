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
decision is made but the refactor has not landed, except in a row marked
"candidate only", which records a need whose crate has not been measured yet.
Crate wave 1 (2026-09-14)
landed three refactors, part of a fourth, and overturned the digamma decision;
each row names the commit that carries it. Early-TOPP wave 2 (2026-09-14)
reopened the Levenberg-Marquardt row after B3-LM's gate and added two "keep"
rows, the quadtree and the overall-score power.

## Replace with a crate

| C++ library and use | Rust today | Crate, pin, features | Status | Measured basis |
|---|---|---|---|---|
| evergreen FFT (`KernelDensityEstimation`) | `src/math/fft.rs` | `rustfft =6.4.1`, always `FftPlannerScalar` | Done, `d396f42` | Upstream KDE expectations pass unchanged. The SIMD planner differs from scalar in 92-99% of components for at most 1.1x on KDE lengths. |
| Eigen `unsupported/NonLinearOptimization` Levenberg-Marquardt (Gauss, Gamma, Gumbel, Gumbel-MLE fitters; `TraceFitter`, later `LevMarqFitter1D`) | `src/math/fitters/levenberg_marquardt.rs` (the Eigen/MINPACK transcription, kept) | `levenberg-marquardt =0.14.0`; `nalgebra =0.33.3`, `default-features = false`, `alloc` + `libm`; pinned in `Cargo.toml`, used only by the ignored candidate test | Reopened: measured, not adopted (B3-LM, `f604b4a`; decision D2's fallback) | Survey: the swapped solver passes all 30 fitter, 9 PEP and 28 unit tests under 1.85, and on the one C++-published Gauss case it is slightly closer to C++ than the port (A 1.5e-12 vs 2.2e-12). B3-LM gate (dax, debug and release identical): an adapter with exact Eigen `maxfev` emulation matches Eigen's status, nfev and njev at 29,003 of 29,004 C2 trace-fit budgets (max_fev 1..500 over the 8 class-test and 50 FeatureFinderCentroided_1 fits, and 4 degenerate fits at 500) and the transcription at all 8,000 distribution-fitter budgets. The miss is `degenerate/flat3_gauss`: status 4 after 24 evaluations against Eigen's 2 after 16, consistent with the crate's exact-zero QR rank against Eigen's thresholded rank (inferred from the paths, not instrumented). The parameter clause fails: at 3,811 of 29,004 trace budgets the candidate's x is neither within 1e-12 of the transcription nor at most as far from C2, 9 of them at budget 500 (FFC_1 Gauss seeds 10, 12 and 20: 2.2e-9 from C2, where the transcription is 6.4e-10); 2,533 of 8,000 distribution budgets lie beyond 1e-12 of the transcription (up to 2.4e-7 for the maximum-likelihood fit), with no C++ sweep to judge them. `minpack-compat` is not closer (3,811 again). MSRV 1.85, pure Rust, `nalgebra` `alloc` + `libm` without `matrixmultiply`; about 1.26x faster per fit in release (1.20 s against 1.51 s for 62 fits x 200 on dax; 1.28x on spock), a gap B3 attributes, without measuring it, to the transcription's per-iteration allocations. The transcription matches C2 in status, nfev and njev at all 29,004 budgets and stays the backend. Lane B3b is root-causing the transcription's own first-step divergence from Eigen. Reopen if a crate reproduces Eigen's thresholded QR rank, its `pnorm / 0.1` delta update and its first-iteration step clamp on every trial until a step is accepted. Whether the two pins move to `[dev-dependencies]` or go with the candidate is open for the lead. Details: [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) section 8. |
| Boost.Math normal quantile (`MultipleTesting` lfdr) | `standard_normal_quantile` in `src/math/multiple_testing.rs` | `statrs =0.18.0`, `default-features = false`; call `erfc_inv` directly | Done, `823c8e9` (merge of `29bfbc4`) | `-SQRT_2 * erfc_inv(2p)` is Boost's own formula, now in Boost's statement order. Measured after the refactor: at most 1.55 ulp from correctly rounded quantiles over the 40 reference points (2.93e-16 relative), where the replaced Acklam+Halley was up to 1.1e-9 relative off; bit-identical to Boost 1.92's double quantile at 79.7% of 97,489 grid points on macOS arm64, against 46.1% before. The survey had measured 4.4e-16 against 1.1e-9 relative to Boost 1.92. lfdr tests unchanged. |
| Boost.Random `mt19937_64` engine (`DecoyGenerator` via `Math::RandomShuffler`; `std::mt19937_64` in `UniqueIdGenerator`) | engine in `src/chemistry/decoy_random.rs` | `rand_mt =6.0.3`, `default-features = false` | Done, `f913cf5` (merge of `1d0e9c8`) | 0 differing words in 52M across 20,012 seeds. The 13 upstream DecoyGenerator strings pass. Boost's `normalize` never affects output. Keep the Boost `uniform_int` bucket mapping and the Fisher-Yates loop: `rand`'s mapping changes every decoy. Refactor check: 0 of 33,310,000 words and 0 of 249,750 shuffles differ from the replaced engine; the engine size stays 2504 bytes. |
| Boost.Regex (`SpectrumLookup`, `SpectrumNativeIDParser`, `MascotXMLFile`, `PepXMLFile`, `MzIdentMLHandler`, `MzTabFile`, `IndexedMzMLDecoder`, `FalseDiscoveryRate`, `EnzymaticDigestion`, `RNaseDigestion`) | Hand-written per-pattern scanners; arbitrary expressions refused with `Unsupported` | `fancy-regex =0.19.2`, behind one Boost-compatible facade | Pending; the facade is in review on `crate/regex-facade` and not merged | The facade needs ASCII bytes mode, dot-matches-newline, `^`/`$` rewritten to Boost's line anchors (`\n`, `\r`, `\f`) and a `match_not_initial_null` retry. With those, 0 mismatches on 240,223 cases over 65 source patterns. The same corpus exposes three port bugs: `mzidentml` `split_annotation` drops multi-loss annotations, `mztab` misses nested brackets, `pepxml` anchors `$` at end of string and parses u64 where C++ uses int32. The plain `regex` crate rejects the lookaround enzyme expressions and duplicate group names. |
| Xerces character/entity reference expansion on read | `expand_reference` in `mascot_xml.rs`; inline copies in `identification_xml.rs` and `mzidentml.rs` | `quick-xml` (existing): `BytesRef::resolve_char_ref`, `escape::resolve_xml_entity` | Partly done, `db19767` (merge of `be2bb15`): `identification_xml.rs` only | 15/20 probes identical. The 5 that differ (`&#X2E;`, `&#+46;`, `&#x+2E;`, `&#0;`, `&#x0;`) are accepted by the hand-rolled code and rejected by quick-xml and libxml2, as the XML 1.0 grammar requires. In `identification_xml.rs` element text now resolves through `resolve_char_ref` and `resolve_xml_entity`; the only accepted-set change there was signed references (`&#+46;`, `&#x+2E;`), because `&#X2E;`, `&#0;` and `&#x0;` were already refused (14,141,319 reference names compared). Attribute values still go through `escape::unescape`, which uses the same character-reference parser. Remaining: `expand_reference` in `mascot_xml.rs` and the inline copy in `mzidentml.rs`. |
| libcurl URL parsing (`NetworkGetRequest`) | host-stripping in `check_http_url` | `ureq::http::Uri` (existing, via ureq) | Done, `db19767` (merge of `be2bb15`) | 9/9 existing cases pass. Over 2M random URLs the new check never rejects a usable host, and it closes 275 empty-host URLs that reached DNS. Keep the scheme gate: `http::Uri` rejects `file://` as malformed. Measured after the refactor on the same 2M URLs: 0 previously refused URLs are accepted, the 275 empty-host URLs fail offline, and 866,415 unparseable URLs are refused by the preflight instead of by ureq, with the same error variant. |
| libc `access(file, R_OK)` (`File::readable`, `File.cpp:506-514`; the TOPPBase input-file checks) | `file::readable` in `src/system/file.rs`: opens a regular file or lists a directory, and answers `false` for anything else without opening it | `rustix` with the `fs` feature: `rustix::fs::access(path, Access::READ_OK)`; no pin chosen | Pending; candidate only, not yet measured against the five tests | Needed for `File::readable` parity. For `-in /dev/null` the executed C++ product SDK exits 4 (INPUT_FILE_EMPTY, "Error: File empty (the file '/dev/null' is empty)"). The port exits 2 with a false "not readable" message, because `cli::input_file_readable` (`src/cli/context.rs:378`) asks `file::readable`, which cannot open a device or FIFO without risking a block. A query that never opens the file reaches the empty-file check, while a mode-000 FIFO given as `-in` must stay exit 2. `Cargo.toml` forbids `unsafe_code` and has no direct `libc` or `rustix` dependency (`libc` is in `Cargo.lock` only transitively), so safe code cannot ask `access` today. Read from the registry source, not measured: `rustix` 1.1.4 declares `rust-version = "1.63"`, and `fs::access` wraps the POSIX `access` call, through `linux-raw-sys` on x86_64 and aarch64 Linux and through `libc` on macOS. Still to measure: MSRV on the Linux build node, the packages added to `Cargo.lock`, and the oracle cases `-in /dev/null` (exit 4) and a mode-000 FIFO as `-in` (exit 2). Evidence: CLI-1 fix round 2 (`86d2733`), integrator request 1 and the verifier's third minor finding. |

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
| extern/Quadtree (Pierre Vigier, MIT; `FeatureOverlapFilter`) | `src/processing/feature_overlap_filter/quadtree.rs` | B9-OVERLAP (`f1c7785`) surveyed seven crates and none reproduces the query order the filter's callbacks observe (16 values per leaf, depth 8, `f32` boxes, strict intersection, a node's values before its children, NW, NE, SW, SE, and the source's split, removal and merge rules): `aabb-quadtree` 0.2.0 (inclusive and epsilon intersection, duplicated items, results sorted by id, unmaintained since 2019), `quadtree-f32` 0.5.0 (four items per node, inclusive overlap), `quadtree` 0.5.0 (points only), `quadtree_rs` 0.1.3 (integer regions), `rectutils` 0.7.0 (static build, entries copied into every leaf), `spart` 0.6.1 (point trees) and `rstar` 0.13 (an R*-tree). The reasons are recorded in the support document and in the provenance `crate_survey`; the package reports no build or measurement of a candidate against the five tests, and the verifier checked the claims for plausibility only, because none of the crates is in the local registry. The port replays the pinned header, executed with and without `NDEBUG`, bit for bit. Reopen if a crate matches that query order exactly. |
| C `powf` through `std::pow(float, float)` (`FeatureFinderAlgorithmPicked` overall seed score) | `overall_score` in `src/analysis/feature_finder_picked/seeds.rs`: `libm::pow` in `f64`, rounded once to `f32` (existing `libm =0.2.16`) | Not a new crate: the platform `powf` is not correctly rounded on macOS (99 of 30,840 executed overall scores one binary32 step off, 12 of 3,084 on FFC_1, checked against 60-digit decimal arithmetic), so bitwise agreement with it cannot hold on every machine. `libm::powf`, the direct binary32 port, was measured and rejected: it differs from the oracle on 226 of 3,084 FFC_1 scores. `libm::pow` in `f64` rounded once gives the correctly rounded value for all 99 and the same bits everywhere, and no seed list changes (B6-FFAP-SEEDS, `80bbdf1`; CPP-272). Following the platform `powf` instead is decision (c) of B6, open for the lead. |
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
- **Evaluation budget.** Measured by B3-LM: `maxfev` emulation around the
  crate is exact (the one mismatch is a path divergence, not a counting
  error); the gate failed on parameter fidelity; see
  [DISTRIBUTION_FITTERS_SUPPORT](DISTRIBUTION_FITTERS_SUPPORT.md) section 8.
- **Feature unification.** `nalgebra/std` enabled anywhere in the graph silently
  brings back `matrixmultiply`'s runtime dispatch. Check `cargo tree -e features`
  when adding dependencies.
