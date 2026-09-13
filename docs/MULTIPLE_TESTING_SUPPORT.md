# MultipleTesting: native header equivalent

[`src/math/multiple_testing.rs`](../src/math/multiple_testing.rs) provides a
native equivalent for the public surface of core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` `MATH/STATISTICS/MultipleTesting.h`
and `source/MATH/STATISTICS/MultipleTesting.cpp`.

Tests: [`tests/multiple_testing.rs`](../tests/multiple_testing.rs). Manifest:
[`tests/data/math_kde_provenance.json`](../tests/data/math_kde_provenance.json).
It builds on [`KERNEL_DENSITY_SUPPORT.md`](KERNEL_DENSITY_SUPPORT.md) and
[`RANK_DATA_SUPPORT.md`](RANK_DATA_SUPPORT.md).

## API mapping

Every public member of the header appears here.

| Source member | Native representation |
| --- | --- |
| `struct Pi0Result` | `Pi0Result` |
| `Pi0Result::pi0` (default `1.0`) | `Pi0Result::pi0`; `Default` carries `1.0` |
| `Pi0Result::pi0_lambda` | `Pi0Result::pi0_lambda` |
| `Pi0Result::lambda_` | `Pi0Result::lambda`; the trailing underscore is a C++ member convention, not part of the name |
| `Pi0Result::pi0_smooth` (default `false`) | `Pi0Result::pi0_smooth` |
| `struct MultipleTesting` | the module; the type holds no state |
| `enum class MultipleTesting::Pi0Method { Smoother, Bootstrap }` | `Pi0Method`, `Default` = `Smoother` |
| `enum class MultipleTesting::LfdrTransform { Probit, Logit }` | `LfdrTransform`, `Default` = `Probit` |
| `static std::string pi0MethodToString(Pi0Method)` | `Pi0Method::as_str() -> &'static str`; the source's `default: return "unknown"` arm cannot be reached from a Rust enum and has no counterpart |
| `static Pi0Method toPi0Method(const std::string&)` | `Pi0Method::parse(&str) -> Result<Self>`, case-insensitive as the source is |
| `static std::string lfdrTransformToString(LfdrTransform)` | `LfdrTransform::as_str()` |
| `static LfdrTransform toLfdrTransform(const std::string&)` | `LfdrTransform::parse(&str) -> Result<Self>` |
| `static std::vector<double> qValue(const std::vector<double>&, double pi0, bool pfdr = false)` | `q_value(&[f64], f64, bool) -> Result<Vec<f64>>` |
| `static Pi0Result pi0Est(const std::vector<double>&, const std::vector<double>& lambda_ = {}, Pi0Method = Smoother, int smooth_df = 3, bool smooth_log_pi0 = false)` | `pi0_est(&[f64], &[f64], Pi0Method, i32, bool, Option<&dyn Pi0Smoother>) -> Result<Pi0Result>`; the extra parameter is the smoothing spline, see differences |
| `static std::vector<double> lfdr(const std::vector<double>&, double pi0, bool trunc = true, bool monotone = true, LfdrTransform = Probit, double adj = 1.5, double eps = 1e-8, std::size_t gridsize = 512, double cut = 3.0)` | `lfdr(&[f64], f64, &LfdrOptions) -> Result<Vec<f64>>`; the six trailing defaults become `LfdrOptions`, whose `Default` is that argument list |
| `static std::vector<double> pNorm(const std::vector<double>& stat, const std::vector<double>& stat0)` | `p_norm(&[f64], &[f64]) -> Result<Vec<f64>>` |
| `template <class T> static std::vector<double> computeModelFDR(const std::vector<T>&)` | `compute_model_fdr(&[f64]) -> Result<Vec<f64>>`; see differences on the template |
| `template <class T> static std::vector<double> pEmp(const std::vector<T>&, const std::vector<T>&)` | `p_emp(&[f64], &[f64]) -> Result<Vec<f64>>` |
| free `qValue(...)` (inline wrapper) | the same `q_value`; a namespace-level alias for a static member is not reproduced |
| free `pi0Est(..., const std::string& pi0_method = "smoother", ...)` | `Pi0Method::parse` plus `pi0_est` |
| free `lfdr(..., const std::string& transf = "probit", ...)` | `LfdrTransform::parse` plus `lfdr` |
| free `pNorm(...)`, `computeModelFDR<T>(...)`, `pEmp<T>(...)` (inline wrappers) | the same functions |
| file-static `argsort_asc` | private `argsort_asc` |
| file-static `percentile` | private `percentile` |
| Not in source but added | `MAX_ITEMS`, `DEFAULT_LFDR_ADJ`, `DEFAULT_LFDR_EPS`, `DEFAULT_SMOOTH_DF`, `LfdrOptions`, `Pi0Smoother` |

## Preserved source conventions

- **`qValue`.** Non-finite entries are dropped, the rest validated to `[0, 1]`,
  and `q = pi0 * m * p / rank(p)` computed with the **`max`** tie rank over the
  `m` finite p-values. A zero rank would give `+inf` rather than a division by
  zero. The monotone sweep caps the largest p-value's q at `1`, replaces a `NaN`
  there with `1.0`, and then walks the p-ordering right to left taking a running
  minimum. The `pfdr` variant divides additionally by `1 - (1-p)^m`.
- **`pi0Est`'s default lambda grid is built by repeated addition** —
  `for (double l = 0.05; l < 1.0 - 1e-12; l += 0.05)` — so it is not exactly
  `0.05 k`: the third entry is `0.15000000000000002`. That is reproduced
  literally, because the smoothing spline is fitted through those abscissae.
  Nineteen entries, `0.05` to `0.95`.
- **A single lambda short-circuits**: `min(#{p >= l} / (m (1-l)), 1)`, with
  `pi0_smooth = false`.
- **The smoother's fallbacks**, all producing `min(min(pi0_lambda), 1)` — the
  *smallest* per-lambda estimate, not a cap — with
  `pi0_smooth = false`: fewer than four lambdas, fewer than two distinct
  lambdas, a spline that does not fit, a `NaN` prediction, and — only when
  `smooth_log_pi0` is off — a non-finite prediction. With `smooth_log_pi0` on, a
  non-finite prediction is *not* a fallback trigger, because `exp(-inf)` is a
  legitimate zero. The successful prediction is clamped to `[0, 1]`.
- **Duplicate lambdas: last value wins.** The source builds a
  `std::map<double, double>` with `xy[lambda] = y`, and `operator[]` assignment
  overwrites. The port inserts into a sorted vector with the same rule.
- **The bootstrap branch** uses the source's own `percentile`, nearest-rank with
  `floor(0.1 (n-1) + 0.5)`, which neither interpolates nor matches
  `Math::quantile`; then minimises
  `W/(m^2 (1-l)^2) (1 - W/m) + (pi0_lambda - minpi0)^2` with a strict `<`, so
  the **first** minimum wins.
- **`pNorm`** computes `1 - 0.5 (1 + erf(z / sqrt 2))` rather than a library
  survival function, and treats a zero-variance null as a point mass: `1.0`
  below the mean, `0.0` at or above it. Non-finite statistics yield `NaN`.
- **`lfdr`'s probit branch clips `p` in place** to `[eps, 1-eps]` *before* the
  transform, so the clipped values are also what the null density is evaluated
  at and what the monotone sort orders by. The logit branch does not clip; `eps`
  goes inside the log-odds instead.
- **`lfdr` divides `pi0 f0` by the density estimate** and yields `+inf` where
  the estimate is not positive, which `truncate` then caps at `1`.
- **`computeModelFDR`'s all-or-nothing `NaN`.** Any `NaN` anywhere returns an
  all-`NaN` vector of the input length and raises nothing. The header documents
  this as differing from `qValue` and warns callers to pre-filter. Reproduced,
  not repaired: a caller reading one position cannot tell a repaired vector from
  a valid one.
- **`pEmp`'s floor.** Everything at or below `1 / |stat0|` is raised to it, with
  `<=` so a value exactly at the floor is rewritten to the same number.
- **`pEmp`'s rank indirection**: `floor(rankdata(-stat, 'average')) - 1`,
  clamped to the last index.

## Native differences

- **`pi0Est` takes the smoothing spline as a parameter.** The source constructs
  `OpenMS::Math::BSplineSmoothingSpline spl(xs, ys, -1.0, smooth_df)` and reads
  `spl.eval(max_lambda)`, guarding on `spl.ok()`. That class is ported, as
  `BSplineSmoothingSpline` in `src/processing/spline/smoothing.rs`, but it lives
  in the `processing` module and `tools/check_module_cycles.py` freezes the
  cross-module edge set: `math` reaches no other top-level module and may not
  start. Duplicating an 1,800-line spline implementation inside `math` would be
  worse than stating the requirement, so `pi0_est` takes
  `Option<&dyn Pi0Smoother>` and `None` behaves exactly as the source's
  `!spl.ok()` fallback. `tests/multiple_testing.rs` implements the trait with
  the ported `BSplineSmoothingSpline` and reproduces the class test's pi0
  literals through it, so the seam carries the source's behaviour rather than
  merely declaring it. Removing the seam needs a restructuring of the module
  graph and is recorded as a deferral.

  **Passing `None` is a different answer, and the less safe one.** In the
  source `!spl.ok()` is a pathology; here it is what a caller who omits the
  smoother gets. Measured on `tests/data/test_lfdr_ref_data.csv`, the 3,170
  PyProphet p-values the class test ships, with the default lambda grid:

  | call | `pi0` |
  | --- | --- |
  | `pi0_est(p, &[], Smoother, 3, false, Some(spline))` | `0.6685639` |
  | `pi0_est(p, &[], Smoother, 3, false, None)` | `0.6403785` |

  The first reproduces the C++ default call's literal `0.6685638`; the second
  is `min(min(pi0_lambda), 1)` and is **lower**. `qValue` computes
  `q = pi0 * m * p / rank(p)` and `lfdr` computes `pi0 * f0 / y`, both linear
  in `pi0`, so the lower estimate makes every q-value and every local FDR
  *smaller* — more hypotheses clear any fixed threshold. Omitting the smoother
  loosens the correction rather than tightening it. `tests/multiple_testing.rs`
  pins both numbers and the direction between them, so the claim is a test and
  not a sentence. `Pi0Result::pi0_smooth` reports which path was taken.
- **`lfdr`'s six trailing parameters become `LfdrOptions`.** Nine positional
  arguments, six of them defaulted, is not a Rust signature.
- **`computeModelFDR` and `pEmp` take `f64` rather than a generic.** Only the
  `double` instantiation is reachable from ported code, and the source's
  `if constexpr (is_floating_point<T>)` `NaN` branch is vacuous for an integral
  `T`. An integral caller widens.
- **`adj` and `eps` are validated**: both must be finite and positive, and `eps`
  below `0.5`. The source checks neither, and a non-positive `eps` would push
  clipped p-values outside `[0, 1]`.
- **The probit quantile comes from a crate and has infinite limits.** Boost's
  `quantile(normal_distribution(0, 1), p)` is reproduced in Boost's own four
  statements over `statrs::function::erf::erfc_inv`, a port of Boost's inverse
  error function, rather than by a hand-written approximation. At `p = 0` and
  `p = 1` it returns the quantile's infinite limits, where Boost raises
  `std::overflow_error`. Its last bits can differ from Boost's, and for rare
  inputs between machines. See "The probit quantile".
- **Ordering uses `f64::total_cmp`** where the source uses `<` inside a
  `std::stable_sort`. The two differ only on `NaN` — which is filtered out
  before every sort in this module — and on `-0.0` versus `0.0`, which `<`
  treats as equal and `total_cmp` orders. Since the comparator feeds a stable
  sort, equal keys keep input order either way, so the only observable
  difference would be a p-value of `-0.0`, which is outside the validated range
  only in sign.
- **Bounded work.** `MAX_ITEMS = 50,000,000` at every entry point.
- **Serial**, as the source is: `MultipleTesting.cpp` carries no `#pragma omp`.

## Source defects found

1. **`lfdr`'s monotone step does not apply the min-rank it claims.** The code
   builds `assigned`/`minrank` and comments "compute min-rank mapping (rankdata
   'min')", but `assigned` is indexed by the *original position* rather than by
   the value. Since `order` is a permutation, every original position is visited
   exactly once and `minrank` is simply the inverse permutation. The mapping
   therefore undoes the sort and nothing else: tied p-values keep the distinct
   values the cumulative maximum left them, where a real min-rank would give
   them all the first one's value. Observable whenever the input has tied
   p-values. The port reproduces the code, not the comment.
2. **`Pi0Result::pi0_lambda` is not filled consistently.** On the single-lambda
   path it holds the clamped estimate; on every other path it holds the raw,
   unclamped ones, which routinely exceed `1`. Reproduced, and stated at the
   field.

Both are reported in the work package's C++ issue list.

## Checked boundaries and evidence

| Boundary | Where |
| --- | --- |
| `MAX_ITEMS = 50,000,000` values | `check_items`, at every entry point |
| p-values in `[0, 1]`, `pi0` in `[0, 1]` | `q_value`, `lfdr`, `pi0_est` |
| lambda in `[0, 1)`, at least one lambda, at least one finite p-value | `pi0_est` |
| `stat0` non-empty and with a finite value | `p_norm` |
| both samples non-empty | `p_emp` |
| `adj` finite and positive, `eps` finite and in `(0, 0.5)` | `lfdr` |
| every division guarded: a zero rank gives `+inf`, a non-positive density gives `+inf`, a zero-variance null becomes a point mass | `q_value`, `lfdr`, `p_norm` |

Evidence:

- **Tier 3** for the class test's literals, section by section, with the source
  line above each Rust test.
- **Tier 3** for the two CSV fixtures the upstream test ships and this port
  copies into `tests/data/`. `test_qvalue_ref_data.csv` is R `qvalue` package
  output over 3,170 p-values and is reproduced by `q_value` for both the FDR and
  the positive-FDR variant to `1e-4`. `test_lfdr_ref_data.csv` is PyProphet
  output for four parameter settings — default, `monotone = false`,
  `logit`, and `eps = 1e-2` — and is reproduced by `lfdr` to `1e-2`, the class
  test's own tolerance. Between them they exercise the whole chain: the FFT, the
  kernel density estimate, the cubic spline and the probit transform.
- **Tier 3** for the pi0 literals `0.697161`, `0.6685638` and `0.6658949`, which
  additionally pin the smoothing-spline seam and the accumulated lambda grid.
- **Tier 4** for the values derived in the test and marked there: the q-values of
  `[0.01, 0.02, 0.03]` shown to be exactly `0.03` from the closed form; the
  `pNorm` tail at the mean being exactly `0.5` by symmetry, and equal statistics
  receiving bit-identical tails; `pEmp`'s floor of `1/|stat0|`; the tied-PEP
  `computeModelFDR` value `0.5/3 = 1/6`; the lambda grid's
  `0.15000000000000002`; monotonicity of the q-values in p.

- **Tier 4** for the probit quantile. The private `standard_normal_quantile` is
  asserted within four machine epsilons relative of correctly rounded quantiles
  at 40 values of `p`. It is also asserted to be exactly antisymmetric at every
  dyadic `p = 2^-k` with `2 <= k <= 53`, and exactly `+0.0` at the median. See
  "The probit quantile" below.

No retained C++ output and no oracle driver exists for this header, so no tier 1
or tier 2 claim is made.

## The probit quantile

`lfdr`'s probit branch calls
`boost::math::quantile(boost::math::normal_distribution<double>(0.0, 1.0), p)`
(`MultipleTesting.cpp:452-453`). Boost evaluates that in four statements,
`boost/math/distributions/normal.hpp:251-254` in Boost 1.92:
`result = erfc_inv(2 * p)`, `result = -result`, `result *= sd * root_two` and
`result += mean`.

`standard_normal_quantile` keeps those four statements with `sd = 1` and
`mean = 0`. `erfc_inv` comes from `statrs =0.18.0` with default features off,
as `statrs::function::erf::erfc_inv`, a port of Boost's inverse error function
approximations. Multiplying by an `sd` of one is exact. Boost's `root_two`
literal rounds to the same `f64` as `std::f64::consts::SQRT_2`. The closing
`+ 0.0` changes only the sign of a zero: `p = 0.5` gives `+0.0`, as Boost does,
and as the replaced code also did.

The crate decision and its survey evidence are recorded in
[`THIRD_PARTY_CRATE_DECISIONS.md`](THIRD_PARTY_CRATE_DECISIONS.md).

### What this replaced

Until this change the port used Acklam's rational approximation with one
Halley step against `libm::erfc`. This document claimed that was accurate to "a
few units in the last place". Measured, it was not:

| Implementation | Worst distance from the correctly rounded quantile, 40 points |
| --- | --- |
| `statrs` `erfc_inv` in Boost's statement order (now) | 1.55 ulp, `2.9e-16` relative |
| Acklam plus one Halley step (before) | `1.1e-9` relative at `p = 1 - 1e-13`; `2.7e-14` at `p = 0.499` |
| `statistics.NormalDist` (AS241), cross-check only | 3.42 ulp |

The measurements ran on x86_64 Linux. The reference quantiles were computed to
110 significant digits by `../oracle/quantile-lane/derive_normal_quantiles.py`,
standard-library Python whose sha256 is in
[`math_kde_provenance.json`](../tests/data/math_kde_provenance.json). The 40
points reach every branch of Boost's `erf_inv` that a binary64 `p` can, both
tails down to `f64::MIN_POSITIVE` and `1 - 1e-15`, and the neighbourhood of the
median. At four epsilons relative, the unit test's bound, the old code fails at
12 of them.

### Output bits changed

These counts compare `lfdr` outputs before and after the change, on x86_64
Linux:

| Input set | Values | Bits changed | Largest change |
| --- | --- | --- | --- |
| Clipped fixture p-values fed to the quantile, `eps = 1e-8` | 3,170 | 2,097 | 934 ulp, `1.6e-13` relative |
| Clipped fixture p-values fed to the quantile, `eps = 1e-2` | 3,170 | 2,237 | 934 ulp, `1.6e-13` relative |
| `lfdr`, default options | 3,170 | 2,237 | `2.1e-15` |
| `lfdr`, `monotone = false` | 3,170 | 2,734 | `3.0e-14` (`8.6e-14` relative) |
| `lfdr`, `eps = 1e-2` | 3,170 | 1,415 | `4.4e-16` |
| `lfdr`, the six-value basic check | 6 | 3 | `2.2e-16` |
| `lfdr`, logit (does not call the quantile) | 3,170 | 0 | none |

Every assertion passes unchanged, including the class test's `1e-2` against
the PyProphet reference.

### Machines

`statrs` has no SIMD and no runtime CPU dispatch; `erf_inv` is scalar Horner
evaluation. It does take `ln` and `sqrt` from the platform's math library, as
the crate's own `.ln()` call sites already do.

Two survey evaluations measured quantiles on aarch64 macOS and x86_64 Linux.
One found 14 of 109,361 inputs differing by one or two ulp, the other 1 of
24,430. None of the differing inputs is a fixture input. The replaced code
differed in none; the survey attributes that to its final Halley step against
the pure-Rust `libm::erfc`. So **`lfdr` results can now differ between machines
in the last place for rare inputs**.

Boost itself is not one bit pattern either. It promotes `double` to 80-bit
`long double` on x86_64 Linux, and its quantile differs between those two
platforms in 31% of the survey's grid.

### Limits

At `p <= 0` and `p >= 1` the function returns the quantile's infinite limits.
Under its default policy, Boost raises `std::overflow_error` at `p = 0` and
`p = 1` instead (`boost/math/special_functions/detail/erf_inv.hpp:363-366`).

`lfdr` never passes `p <= 0`, because it clips to `eps` and validates `eps` as
positive. It passes `p = 1` only when `eps` is at most `2^-54` (about
`5.6e-17`), where `1 - eps` rounds to one. There the source would throw and the
port returns `+inf`, exactly as it did before this change.
