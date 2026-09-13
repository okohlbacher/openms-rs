# MATH/STATISTICS distribution fitters

Port of the four `MATH/STATISTICS` fitter headers at
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`:

| C++ header | C++ implementation | Rust |
|---|---|---|
| `MATH/STATISTICS/GaussFitter.h` | `GaussFitter.cpp` | `src/math/fitters/gauss.rs` |
| `MATH/STATISTICS/GammaDistributionFitter.h` | `GammaDistributionFitter.cpp` | `src/math/fitters/gamma.rs` |
| `MATH/STATISTICS/GumbelDistributionFitter.h` | `GumbelDistributionFitter.cpp` | `src/math/fitters/gumbel.rs` |
| `MATH/STATISTICS/GumbelMaxLikelihoodFitter.h` | `GumbelMaxLikelihoodFitter.cpp` | `src/math/fitters/gumbel_max_likelihood.rs` |

All four share `src/math/fitters/levenberg_marquardt.rs`, which has no C++
header of its own: it reproduces the `Eigen::LevenbergMarquardt` the four
`.cpp` files call. `src/math/mod.rs` and `src/math/fitters/mod.rs` are the new
domain roots. Tests are `tests/math_distribution_fitters.rs` plus the
`#[cfg(test)]` modules inside each source file; the manifest is
`tests/data/distribution_fitters_provenance.json`.

The four headers and their implementations are 928 physical lines (385 of
header, 543 of implementation) and their class tests carry 25 `START_SECTION`
blocks. All 25 are mapped; see **Class-test sections** below.

---

## 1. Why the solver is ported too

The C++ does not implement an optimizer. `GaussFitter`,
`GammaDistributionFitter` and `GumbelDistributionFitter` each build a functor
with an analytic Jacobian and hand it to `Eigen::LevenbergMarquardt<Functor>`;
`GumbelMaxLikelihoodFitter` wraps its functor in `Eigen::NumericalDiff` first.
Eigen's class is a transcription of MINPACK `lmder`, and the crate may not take
a dependency, so it is reproduced in
`src/math/fitters/levenberg_marquardt.rs`:

* column-pivoted Householder QR with LAPACK's norm-downdating rule and its
  `sqrt(eps)` recomputation threshold, and Eigen's rank threshold
  `|maxpivot| * eps * min(rows, cols)`;
* `lmpar` in Eigen's `lmpar2` form, capped at ten iterations, with `qrsolv`
  eliminating the scaled diagonal by Givens rotations;
* MINPACK's termination ladder, in order: `RelativeErrorAndReductionTooSmall`,
  `RelativeReductionTooSmall`, `RelativeErrorTooSmall`,
  `TooManyFunctionEvaluation`, `FtolTooSmall`, `XtolTooSmall`, `GtolTooSmall`;
* Eigen's defaults, which none of the four fitters overrides: `factor = 100`,
  `maxfev = 400`, `ftol = xtol = sqrt(f64::EPSILON)`, `gtol = 0`;
* **all three** of the norms this path uses, each where Eigen uses it. An
  earlier revision of this document said "both of Eigen's norms"; that was
  wrong, and the code was not - there are three:

  | Norm | Eigen expression | Where |
  |---|---|---|
  | `stable_norm` | `stableNorm()` - one scaling pass by the largest magnitude, then the scaled sum of squares | `fnorm`, `fnorm1`, `pnorm`, `xnorm`, the predicted-reduction norm, and `gnorm` inside `lmpar` |
  | `blue_norm` | `blueNorm()` - Blue's three-bin algorithm | `fjac.colwise().blueNorm()`, the column norms the driver turns into `diag`; and `dxnorm` and the two correction norms inside `lmpar` |
  | `plain_norm` | `MatrixBase::norm()` - `sqrt` of a plain sum of squares, no scaling pass | the pivot column norms inside `ColPivHouseholderQR`, both the initial pass and the direct recomputation the LAPACK downdating rule falls back to |

  The distinction is load-bearing: substituting any one for another changes the
  pivot order or the trust-region scaling and therefore the answer.

This is what makes the published parameters reachable. A generic
Levenberg-Marquardt with a different stopping rule stops at a different point:
Eigen's `ftol` of `1.49e-8` is a loose relative reduction, so the reported
parameters are a property of *that* stopping rule and not only of the residual
surface.

### The one solver deviation found by audit, and what it cost

An audit of this port against `unsupported/Eigen/src/NonLinearOptimization/`
found one place where the transcription had silently normalised an asymmetry in
the reference, and the earlier revision of this document did not mention it.

`lmpar` forms the quantity `P^-1 (diag .* diag .* x) / dxnorm` twice, and Eigen
spells the two occurrences with **different associations**:

| Site | Eigen, `lmpar.h` | Meaning |
|---|---|---|
| the `parl` lower bound, `rank == n` branch | line 198: `P^-1 * diag.cwiseProduct(wa2) / dxnorm` | `(diag[i] * wa2[i]) / dxnorm` - multiply, then divide |
| the Newton correction inside the loop | line 241: `P^-1 * diag.cwiseProduct(wa2 / dxnorm)` | `diag[i] * (wa2[i] / dxnorm)` - divide, then multiply |

The port had taken the divide-first spelling at both sites, so it agreed with
Eigen's line 241 but not with its line 198. Both lines carry those numbers in
Eigen 3.4.0 - the version installed in the OpenMS build environment - and in
5.0.1, and the asymmetry is identical in both, so it is not an artefact of one
Eigen release. (Whether MINPACK's own Fortran is symmetric here was not
checked: no MINPACK source is available in this tree, and the reference the
port must match is Eigen's C++, which is what the four `.cpp` files call.)

**Which operation, and how large.** One `f64` multiply-divide pair per
parameter, once per `lmpar` call, in the branch that computes the lower bound
`parl`. `(a*b)/c` and `a*(b/c)` differ by at most one unit in the last place of
the result, and only when the two roundings fall on opposite sides; for the
2- and 3-parameter problems here that is a single rounding per parameter per
call, never an accumulation within one call.

**Whether it compounds.** Yes, through the iteration, though not without bound.
`parl` clamps the Levenberg parameter `par` from below, `par` scales the
diagonal handed to `qrsolv`, and the resulting step moves `x`; the next
iteration re-evaluates the residual and the Jacobian at the moved `x`, so a
one-ulp difference re-enters as an O(1 ulp) relative difference in every
subsequent quantity. Measured on the gate node by correcting the association
and re-running the class-test cases:

| Case | Shift caused by the correction |
|---|---|
| `GaussFitter::fit` case 1, `A` | 153 ulp, `3.3e-14` relative |
| `GaussFitter::fit` case 1, `x0` | 21 ulp, `3.9e-15` relative |
| `GaussFitter::fit` case 1, `sigma` | 810 ulp, `1.7e-13` relative |
| `GaussFitter::fit` case 2 | bit-identical either way |
| Gamma, both Gumbel cases, both maximum-likelihood cases | bit-identical either way |

So a single ulp at the start of the long case grows to roughly 800 ulp by the
end of it, and vanishes entirely on the cases that converge in few steps. It is
not visible at the tolerance the class tests assert - `1e-5`, nearly eight
orders of magnitude away - nor at the `1e-9` and `1e-11` this port's own fit
tests assert. Reported for honesty about the margin, not because any assertion
turned on it.

**Status: corrected.** Each site now carries the association of the Eigen line
it transcribes, with that line cited in a comment, and
`the_two_lmpar_scalings_are_not_the_same_association` pins a triple at which
the two spellings disagree, so collapsing them again fails a test. The fit
deviations against the published C++ parameters in §5 were re-measured after
the correction and are marginally smaller for two of the three case-1
parameters and marginally larger for the third - the correction is inside the
`4e-11` noise floor of that case and does not explain it.

**A second, smaller omission found in the same pass, also corrected.**
`stable_norm` took the plain reciprocal `1.0 / max_coeff`, where Eigen's
`stable_norm_kernel` guards both ends of the range first: if `1 / maxCoeff`
overflows - which it does for a subnormal largest coefficient - Eigen uses
`(scale, invScale) = (1 / highest(), highest())`, and if `maxCoeff` is infinite
it uses `(maxCoeff, 1)`. Without the guard the subnormal case returns infinity
instead of a subnormal. Unreachable from the four fitters, which validate their
inputs finite before the solver sees them, but it was a divergence from the
reference and is now transcribed and tested
(`stable_norm_guards_the_reciprocal_at_both_ends_of_the_range`). No published
value changes.

**One accumulation-order difference remains and is not correctable here.**
Eigen's `squaredNorm()` and `dot()` are vectorized reductions: they accumulate
into several SIMD lanes and combine at the end, so their summation order
depends on the compiler, the target and the alignment of the data. This port
accumulates sequentially. Every `squaredNorm()` and `dot()` on this path is
affected: the QR's pivot column norms, the tail norm inside `makeHouseholder`,
the scaled inner sum inside `stableNorm`, and the gradient dot product in
`minimizeOneStep`. `blueNorm` is *not* affected - Eigen's `blueNorm_impl` is a
scalar loop over the coefficients, which is what this port reproduces. The
difference is at most one unit in the last place per reduction, it is not
reproducible across builds of the C++ itself, and it is the leading candidate
for the residual `4e-11` in the first `GaussFitter` case. Matching it would mean
guessing a particular Eigen build's vector width, which would be a worse kind of
infidelity.

---

## 2. API mapping

Every public member of the four headers appears here.

### `GaussFitter.h`

| C++ member | Rust |
|---|---|
| `struct GaussFitter::GaussFitResult` | `gauss::GaussFitResult` |
| `GaussFitResult()` | `GaussFitResult::default()`, `(-1, -1, -1)` |
| `GaussFitResult(double a, double x, double s)` | `GaussFitResult::new(a, x0, sigma)` |
| `double A` | field `a` |
| `double x0` | field `x0` |
| `double sigma` | field `sigma` |
| `double eval(double x) const` | `GaussFitResult::eval(&self, f64) -> Result<f64>` |
| `double log_eval_no_normalize(double x) const` | `GaussFitResult::log_eval_no_normalize(&self, f64) -> Result<f64>` |
| private `double halflogtwopi` | module constant `HALF_LOG_TWO_PI`; see §4 |
| `GaussFitter()` | `GaussFitter::new()` / `GaussFitter::default()` |
| `virtual ~GaussFitter()` | not ported: no owned resources, no `Drop` |
| `void setInitialParameters(const GaussFitResult&)` | `GaussFitter::set_initial_parameters` |
| `GaussFitResult fit(std::vector<DPosition<2>>&) const` | `GaussFitter::fit(&self, &[(f64, f64)]) -> Result<GaussFitResult>` |
| `static std::vector<double> eval(const std::vector<double>&, const GaussFitResult&)` | `GaussFitter::eval(&[f64], &GaussFitResult) -> Result<Vec<f64>>` |
| protected `GaussFitResult init_param_` | private field, readable through `initial_parameters()` |
| private `GaussFitter(const GaussFitter&)` (declared, undefined) | not ported: the port derives `Copy`; see §4 |
| private `GaussFitter& operator=(const GaussFitter&)` (declared, undefined) | not ported: the port derives `Copy`; see §4 |

### `GammaDistributionFitter.h`

| C++ member | Rust |
|---|---|
| `struct GammaDistributionFitter::GammaDistributionFitResult` | `gamma::GammaDistributionFitResult` |
| `GammaDistributionFitResult(double bIn, double pIn)` | `GammaDistributionFitResult::new(b, p)`; no `Default`, as in C++ |
| `double b` | field `b` |
| `double p` | field `p` |
| `GammaDistributionFitter()` | `GammaDistributionFitter::new()` / `default()` |
| `virtual ~GammaDistributionFitter()` | not ported: no owned resources, no `Drop` |
| `void setInitialParameters(const GammaDistributionFitResult&)` | `GammaDistributionFitter::set_initial_parameters` |
| `GammaDistributionFitResult fit(const std::vector<DPosition<2>>&) const` | `GammaDistributionFitter::fit(&self, &[(f64, f64)]) -> Result<GammaDistributionFitResult>` |
| protected `GammaDistributionFitResult init_param_` | private field, readable through `initial_parameters()` |
| private copy constructor and `operator=` (declared, undefined) | not ported: the port derives `Copy` |
| - | `GammaDistributionFitResult::eval` is native: the density lives only inside the C++ functor |

### `GumbelDistributionFitter.h`

| C++ member | Rust |
|---|---|
| `struct GumbelDistributionFitter::GumbelDistributionFitResult` | `gumbel::GumbelDistributionFitResult` |
| `GumbelDistributionFitResult(double a = 1.0, double b = 2.0)` | `GumbelDistributionFitResult::new(a, b)` and `Default` giving `(1.0, 2.0)` |
| `double a` | field `a` |
| `double b` | field `b` |
| `double eval(double) const` | **not ported: declared in the header and defined nowhere in the SDK.** `gumbel::GumbelDistributionFitResult::eval` exposes the residual model the functor uses, which is what the declaration was evidently meant to be, and is labelled native at the item |
| `double log_eval_no_normalize(double) const` | `GumbelDistributionFitResult::log_eval_no_normalize(&self, f64) -> Result<f64>` |
| `GumbelDistributionFitter()` | `GumbelDistributionFitter::new()` / `default()` |
| `virtual ~GumbelDistributionFitter()` | not ported: no owned resources, no `Drop` |
| `void setInitialParameters(const GumbelDistributionFitResult&)` | `GumbelDistributionFitter::set_initial_parameters` |
| `GumbelDistributionFitResult fit(std::vector<DPosition<2>>&) const` | `GumbelDistributionFitter::fit(&self, &[(f64, f64)]) -> Result<GumbelDistributionFitResult>` |
| `GumbelDistributionFitResult fitWeighted(const std::vector<double>&, const std::vector<double>&)` | **not ported: declared in the header and defined nowhere in the SDK.** The class test reaches a weighted fit through `GumbelMaxLikelihoodFitter` instead |
| protected `GumbelDistributionFitResult init_param_` | private field, readable through `initial_parameters()` |
| private copy constructor and `operator=` (declared, undefined) | not ported: the port derives `Copy` |

### `GumbelMaxLikelihoodFitter.h`

| C++ member | Rust |
|---|---|
| `struct GumbelMaxLikelihoodFitter::GumbelDistributionFitResult` | `gumbel_max_likelihood::GumbelDistributionFitResult`, a type distinct from the one above |
| `GumbelDistributionFitResult(double a, double b)` | `GumbelDistributionFitResult::new(a, b)`; no `Default`, as in C++ |
| `double a` | field `a` |
| `double b` | field `b` |
| `double log_eval_no_normalize(double) const` | `GumbelDistributionFitResult::log_eval_no_normalize(&self, f64) -> Result<f64>` |
| `GumbelMaxLikelihoodFitter()` | `GumbelMaxLikelihoodFitter::new()` / `default()` |
| `GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)` | `GumbelMaxLikelihoodFitter::with_initial_parameters` |
| `virtual ~GumbelMaxLikelihoodFitter()` | not ported: no owned resources, no `Drop` |
| `void setInitialParameters(const GumbelDistributionFitResult&)` | `GumbelMaxLikelihoodFitter::set_initial_parameters` |
| `GumbelDistributionFitResult fitWeighted(const std::vector<double>&, const std::vector<double>&)` | `GumbelMaxLikelihoodFitter::fit_weighted(&mut self, &[f64], &[f64]) -> Result<GumbelDistributionFitResult>` |
| protected `GumbelDistributionFitResult init_param_` | private field, readable through `initial_parameters()`; overwritten by a successful fit, as in C++ |
| private copy constructor and `operator=` (declared, undefined) | not ported: the port derives `Copy` |

### Native additions

`initial_parameters()` on all four fitters (the C++ member is protected with no
accessor), `gamma::GammaDistributionFitResult::eval`,
`gumbel::GumbelDistributionFitResult::eval`, and the whole of
`levenberg_marquardt`: `MAX_POINTS`, `MAX_BYTES`, `preflight_points`,
`LmParameters`, `LmStatus`, `LmStatus::code`, `DenseMatrix` and its accessors,
`stable_norm`, `blue_norm`, `minimize` and `numerical_jacobian`.

---

## 3. Preserved source conventions

* **The initial guesses, verbatim.** `GaussFitter` starts at
  `(A, x0, sigma) = (0.06, 3.0, 0.5)` (`GaussFitter.cpp:23`),
  `GammaDistributionFitter` at `(b, p) = (1.0, 5.0)`
  (`GammaDistributionFitter.cpp:27`), `GumbelDistributionFitter` and
  `GumbelMaxLikelihoodFitter` at `(a, b) = (0.25, 0.1)`
  (`GumbelDistributionFitter.cpp:32`, `GumbelMaxLikelihoodFitter.cpp:112`).
  None of these is a neutral starting point and the surfaces are not convex.
  An earlier revision asserted, without checking, that a different start "would
  not reach" the published answer of the first `GaussFitter` case. The measured
  situation on that case is more specific and more useful:

  | Start `(A, x0, sigma)` | Outcome |
  |---|---|
  | `(0.06, 3.0, 0.5)`, the source's | `(1.0189827566, 0.3006128709, 0.1363163309)` |
  | `(1.0, 0.3, 0.2)` | same basin, `1.0e-5` relative away in `A` |
  | `(1.0, 1.0, 1.0)` | same basin, `4.6e-6` relative away in `A` |
  | `(0.5, -1.0, 2.0)` | same basin, `1.0e-5` relative away in `A` |
  | `(0.06, 10.0, 0.5)` | a *different* minimum: `A = 1328.4`, `x0 = -582.9`, `sigma = 127.5` |
  | `(0.06, 0.0, 0.5)` | no answer: `TooManyFunctionEvaluation`, which `GaussFitter` turns into `Exception::UnableToFit` |

  So a far-away start can land on another minimum or fail outright, and even a
  start in the right basin lands `1e-5` to `1e-6` relative from the published
  parameters - at least five orders of magnitude further than any arithmetic
  difference discussed in this document. The published numbers belong to the
  source's guess and to no other.
* **The result structs' defaults.** `GaussFitResult` defaults to
  `(-1, -1, -1)` (`GaussFitter.h:42-43`), which is not a usable model.
  `GumbelDistributionFitter`'s result defaults to `(1.0, 2.0)`
  (`GumbelDistributionFitter.h:40`), which is *not* the fitter's starting guess
  of `(0.25, 0.1)` - the constructor assigns over the default-constructed
  member. Both are reproduced, but only one of them is *asserted* by a class
  test, and an earlier revision of this document claimed both were. The correct
  statement: `GumbelDistributionFitter_test.cpp:130-134` has a
  `START_SECTION((GumbelDistributionFitResult()))` that asserts `a == 1.0` and
  `b == 2.0`. `GaussFitter_test.cpp` has no section for the default constructor
  at all; the only place `(-1, -1, -1)` appears there is line 127, where the
  `setInitialParameters` section builds a result through the *three-argument*
  constructor and the section ends in `NOT_TESTABLE`. The port's default is
  therefore evidence-tier "header read", not "class test", and
  `the_default_result_is_the_sources_invalid_marker` is a native test, not a
  transcription.
* **The residual expressions, in the source's arithmetic order.** The Gaussian
  residual reuses `sig2 = 2 * sig * sig` and writes
  `A * exp(-(x - x0) * (x - x0) / sig2) - y` (`GaussFitter.cpp:56`). The Gamma
  residual is `pow(b, p) / tgamma(p) * pow(x, p - 1) * exp(-b * x) - y`
  (`GammaDistributionFitter.cpp:65`). The Gumbel residual forms
  `z = exp((a - x) / b)` first and then `(z * exp(-z)) / b - y`
  (`GumbelDistributionFitter.cpp:65-66`). Reassociating any of these changes the
  last bits and, through the iteration, the reported parameters.
* **The analytic Jacobians, including their redundant sub-expressions.** For
  example `GaussFitter.cpp:77` writes the centre derivative as
  `A * j0 * (-(-2 * x + 2.0 * x0) / sig2)`, which is a double negation of
  `A * j0 * 2 * (x - x0) / sig2`; the port keeps the source's spelling.
* **The customized Gamma density.** `GammaDistributionFitter.h:28-30` documents
  that the fitted function is zero whenever `b <= 0` or `p <= 0`, so that an
  unconstrained optimizer can be used. Both the residual branch
  (`GammaDistributionFitter.cpp:60-74`) and the all-zero Jacobian branch
  (`GammaDistributionFitter.cpp:89-113`) are reproduced.
* **The parameter transforms.** `GaussFitter::fit` reports `|sigma|`
  (`GaussFitter.cpp:107`, with the source's comment that the absolute value is
  the correct solution); `GumbelMaxLikelihoodFitter` uses `|b|` inside the
  objective and reports `|b|` (`GumbelMaxLikelihoodFitter.cpp:56`, `:101`).
  `GumbelDistributionFitter::fit` applies no such transform and its result can
  carry a negative scale; the port does not add one.
* **The status rejection rules.** `GaussFitter` throws for
  `ImproperInputParameters` (0) *and* `TooManyFunctionEvaluation` (5)
  (`GaussFitter.cpp:101-105`); the other three throw for `status <= 0`
  (`GammaDistributionFitter.cpp:132`, `GumbelDistributionFitter.cpp:109`,
  `GumbelMaxLikelihoodFitter.cpp:87`), and after a completed `minimize` the only
  reachable member of that range is `ImproperInputParameters`.
  `LmStatus::code()` exposes the numeric value so the distinction stays legible.
* **The maximum-likelihood objective and its accumulation order.**
  `GumbelMaxLikelihoodFitter.cpp:53-68` accumulates
  `w_i * (-ln sigma - z_i - exp(-z_i))` over the samples in index order into a
  running sum and negates once at the end. The port sums in the same order.
* **Two residuals for a one-dimensional objective.** The same functor declares
  `Functor<double>(2, 2)` (`GumbelMaxLikelihoodFitter.cpp:48`) and leaves
  `fvec(1)` permanently zero, because Eigen rejects a problem with fewer
  residuals than parameters. `RESIDUALS = 2` in the port for the same reason.
* **The forward-difference Jacobian.** The same fitter wraps its functor in
  `Eigen::NumericalDiff` (`GumbelMaxLikelihoodFitter.cpp:80`), so the port uses
  `numerical_jacobian` with Eigen's step, `sqrt(f64::EPSILON) * |x_j|` falling
  back to `sqrt(f64::EPSILON)`, and charges its `n + 1` evaluations against
  `max_fev` exactly as Eigen does. An analytic Jacobian would stop the iteration
  somewhere else.
* **The write-back.** A successful `fitWeighted` overwrites the fitter's own
  start parameters with the result (`GumbelMaxLikelihoodFitter.cpp:98-99`), so a
  second call continues from the first. `fit_weighted` takes `&mut self` for
  that reason, and a rejected call leaves them untouched.
* **Boost's expression order, and the constant Boost actually uses.**
  `GaussFitResult::eval` is `pdf(x) * (A / pdf(x0))`, with the source's comment
  that "simply multiplying the CDF with A is wrong"
  (`GaussFitter.cpp:124`, `:135`). The density is Boost's
  `normal_distribution` pdf in Boost's order: form the deviation, negate and
  square it in place, divide by `2 * sd * sd`, exponentiate, divide by
  `sd * sqrt(2 * pi)`.

  That last divisor was got wrong, and the reason given for it was wrong in the
  specific way this project keeps hitting. An earlier revision of this document
  said the pdf divides by "Boost's `root_two_pi` literal", and the port
  therefore carried `2.506_628_274_631_000_7`, the `f64` that literal rounds to.
  Boost's `normal.hpp` does not use `root_two_pi`: the last line of `pdf` reads

      result /= sd * sqrt(2 * constants::pi<RealType>());

  so the divisor is the square root of the *rounded* `2 * pi`, evaluated at run
  time. For `double` that is `2.506_628_274_631_000_2`, one unit in the last
  place *below* the `root_two_pi` literal. The port now uses that value, named
  `SQRT_TWO_PI`, pinned by a unit test to `(2.0 * PI).sqrt()` and to being one
  ulp below the literal. The correction is worth what a correction of this size
  is ever worth: `GaussFitter::eval` now reproduces all seven published C++
  intensities **bit for bit**, where before it matched four of seven and missed
  the other three by up to `1.8e-16` relative. That `1.8e-16` row in §5 was
  never a `libm` difference; it was this constant.

---

## 4. Native differences

Each is documented at the Rust item as well.

* **`Result` instead of exceptions.** Every `Exception::UnableToFit` becomes
  `Error::InvalidValue` carrying the source's message text, prefixed with the
  source's `UnableToFit-<Class>` tag. The crate's error enum is fixed and has no
  fitting-specific variant.
* **`&[(f64, f64)]` instead of `std::vector<DPosition<2>>`.** The crate has a
  `DPosition` type, in `src/data_structures`, but `src/math` may not depend on
  it: `tools/check_module_cycles.py` freezes the module-pair graph and
  `math -> data_structures` is not in it. A pair of `f64` carries exactly the
  same information for these four call sites. `GaussFitter::fit` and
  `GumbelDistributionFitter::fit` take their points by non-const reference in
  C++ and modify nothing; the port takes a shared slice.
* **Degenerate input is refused rather than answered with NaN.** The source
  checks nothing before constructing the functor. The port rejects, before
  allocating anything: fewer points than parameters (which Eigen would report as
  `ImproperInputParameters` anyway, so the outcome matches, only earlier and
  with a clearer message); non-finite coordinates; negative abscissae for the
  Gamma density, whose `pow(x, p - 1)` is NaN for a negative base; a non-finite
  initial guess; and a zero initial `sigma` or `b`, which would make the first
  residual divide by zero. It also refuses to *return* non-finite parameters.
* **Mismatched weight and sample lengths are an error.**
  `GumbelMaxLikelihoodFitter.cpp:58-63` advances a weight iterator in lockstep
  with the sample iterator and never compares the two sizes, so a short weight
  vector is read past its end. `fit_weighted` requires equal lengths.
* **`eval` and `log_eval_no_normalize` validate.** The source passes `sigma` to
  `boost::math::normal_distribution`, whose default policy throws
  `std::domain_error` for a non-positive or non-finite scale, so `eval` already
  had an error path; `log_eval_no_normalize` did not, and returns NaN for
  `sigma <= 0`. Both return `Result` here, and both also reject a non-finite
  evaluation point.
* **The fitters are `Copy`.** All four C++ classes declare a private copy
  constructor and assignment operator and define neither, which makes them
  non-copyable. The Rust types hold two or three `f64` and nothing else; there
  is no invariant a copy could break, so they derive `Clone, Copy`. The
  class-test sections for those two members are `NOT_TESTABLE` in C++ and are
  mapped to tests that assert a copy really carries the initial guess.
* **`halflogtwopi` is a constant.** In C++ it is a non-static data member with an
  initializer (`GaussFitter.h:76`), so every `GaussFitResult` carries its own
  copy of `0.5 * log(2 * PI)` and the struct is 32 bytes rather than 24. The
  Rust struct is 24 bytes and the value is a module constant. The numeric value
  is identical and a unit test asserts it equals `0.5 * (2.0 * PI).ln()`.
* **`digamma` is not Boost's.** Boost uses rational minimax approximations; the
  port uses the recurrence `psi(x) = psi(x + 1) - 1/x` up to `x >= 10` followed
  by the asymptotic series through `B14`. Measured against the closed forms
  `psi(1) = -gamma` and `psi(1/2) = -gamma - 2 ln 2` and against the recurrence
  at eleven integer points, the implementation agrees to better than `1e-14`
  absolute.

  An earlier revision justified this with "the value enters only the Jacobian,
  which steers the search and does not define the optimum". That is true of the
  minimizer of the sum of squares and false of what `fit` returns, and it
  contradicts the argument §1 makes at length: in Levenberg-Marquardt the
  Jacobian sets `diag`, the trust-region radius, the gradient test and every
  termination test, so the *reported* parameters are a property of the Jacobian
  as well as of the residual. A perturbed digamma can and in general does move
  them.

  The reason the substitution is acceptable is the size of the slack, not an
  absence of effect. `GammaDistributionFitter` is the only one of the four that
  calls digamma, its class test asserts only the parameters the data were
  generated from - `b = 7.25`, `p = 3.11` - at `0.01` absolute, and the port
  lands `2.7e-3` and `6.9e-3` away, with three orders of magnitude between the
  `1e-14` digamma error and that margin. No C++-produced Gamma parameter is
  published anywhere in the SDK, so a tighter claim is not available and is not
  made: this item is an accepted, bounded divergence, not a proven identity.
* **`pow(v, 2.0)` is written `v * v`.** The source spells squares with
  `std::pow` in three places, not the two an earlier revision listed:
  `GaussFitter.cpp:144` (`pow((x - x0) / sigma, 2.0)`),
  `GammaDistributionFitter.cpp:100` (`pow(tgamma(p), 2)`) and
  `GumbelDistributionFitter.cpp:81-85` (`pow(z, 2)` and three `pow(b, 2)`).
  `v * v` is the correctly rounded square by construction; `std::pow(v, 2.0)`
  equals it only if the implementation is correctly rounded at that argument,
  which neither IEEE-754 nor the C standard requires of `pow` in general.
  glibc and Apple's libm both special-case small integral exponents and return
  the rounded product, so on the platforms this crate is gated on the two agree;
  an exotic libm could differ by one unit in the last place. Stated as the
  bounded assumption it is rather than as an identity.
* **`-1.0 * v` is written `-v`.** IEEE negation is exact, so these are the same
  value; the change is only to satisfy `clippy::neg_multiply`.
* **The solver's matrix type is `DenseMatrix`, not `Matrix`.** It is Eigen's
  `MatrixXd` for this one solver and not a port of `DATASTRUCTURES/Matrix.h`;
  the coverage ledger maps candidate Rust types to headers by name, and the
  obvious name would have claimed a header this group does not touch.
* **`lmpar2`'s `s` workspace is `n`-by-`n`.** Eigen copies the whole
  `m`-by-`n` QR factor although `qrsolv` reads only the leading `n`-by-`n` block
  and its diagonal. The port copies the block. The arithmetic is identical; the
  allocation is `m/n` times smaller, which matters because `m` is the number of
  observations.
* **Bounded work.** The C++ has no ceiling. `preflight_points` rejects more than
  `MAX_POINTS = 1_000_000` observations, and a dense Jacobian above
  `MAX_BYTES = 64 MiB`, before anything is allocated.
* **Serial.** Neither the four `.cpp` files nor Eigen's non-linear optimization
  module carries `#pragma omp`, so there is no OpenMP gap to record: the source
  is serial here and so is the port.

  Nothing here is parallelised and nothing here should be. A single fit is a
  sequential dependency chain - iteration `k + 1` cannot start until `k`'s step
  is accepted - and the residual and Jacobian loops are a few dozen to a few
  thousand elements, far below the point where a `rayon` split would pay for
  itself. Recorded as a candidate for a caller, not for this module: fitting
  *many independent* data sets is embarrassingly parallel, because `fit` takes
  `&self` and `GaussFitter`, `GammaDistributionFitter` and
  `GumbelDistributionFitter` are `Copy` with no interior mutability, so a
  `par_iter().map(|d| fitter.fit(d))` over a slice of data sets is already
  sound and, collected in order, bit-identical to the serial loop: each fit
  reduces only over its own data, so there is no cross-item float reduction
  whose order could change. `GumbelMaxLikelihoodFitter::fit_weighted` is the
  exception: it takes
  `&mut self` because the source writes the result back into its own start
  parameters, so a caller wanting that in parallel needs one fitter per task,
  which changes nothing numerically as long as each starts from the same guess.

---

## 5. Checked boundaries and evidence

### Evidence tier

**Tier 3, source review**, for every parameter value: the expected numbers are
transcribed from the four `_test.cpp` files, which are the only oracle available
without building C++. `tests/math_distribution_fitters.rs` marks each
transcription with its source line range.

**Tier 4, independently derived**, for these, which do not depend on any
transcribed literal:

* `mle_result_is_a_local_minimum_of_the_weighted_log_likelihood` recomputes the
  source's weighted negative log-likelihood at the fitted parameters and at
  eight neighbours (`±1e-3` and `±1e-2` in each coordinate) and asserts the
  fitted point is lower than all of them, and that the objective is strictly
  positive there. That is the condition under which minimizing the square of a
  scalar residual is minimizing the residual, so it establishes that the source's
  unusual formulation really does return the maximum-likelihood estimate on this
  data - a claim the transcribed `a ~ 2, b ~ 0.6` at `0.1` tolerance is far too
  loose to support on its own.
* `a_linear_least_squares_problem_is_solved_exactly` drives the solver onto a
  problem with a zero-residual exact solution (`y = 1 + 2x` through three
  points) and asserts it to `1e-10`.
* `blue_norm_constants_are_the_powers_of_two_eigen_derives` asserts each of the
  four Blue's-algorithm constants equals `2.0f64.powi(e)` for the exponent
  Eigen's integer arithmetic produces. Two of the four literals were wrong when
  first transcribed and this test is what caught them.
* `stable_norm_survives_magnitudes_that_overflow_a_naive_sum` checks both norms
  at `1e200` and `1e-200`, where a naive sum of squares overflows or underflows.
* `the_density_divisor_is_the_square_root_boost_evaluates` asserts that the
  density's divisor is exactly `(2.0 * PI).sqrt()`, the expression Boost's
  `normal.hpp` evaluates, and that it is one unit in the last place below
  Boost's unused `root_two_pi` literal. This is what turned `GaussFitter::eval`
  from four-of-seven exact into seven-of-seven; see §3.
* `the_two_lmpar_scalings_are_not_the_same_association` asserts that
  `(a*b)/c != a*(b/c)` at a concrete triple, pinning the asymmetry Eigen's
  `lmpar.h` has between its lines 198 and 241; see §1.
* `the_published_gauss_case_belongs_to_the_sources_initial_guess` re-fits the
  first `GaussFitter` case from five other starting points and asserts that one
  reaches a different minimum, one cannot be fitted at all, and the three that
  reach the right basin still land between `1e-7` and `1e-3` relative from the
  published amplitude. This is the measurement behind the table in §3 and the
  reason the initial guesses and the stopping rule are reproduced verbatim.
* `digamma_reproduces_its_closed_form_values` (closed forms and the recurrence).
* `gauss_log_eval_is_the_log_of_the_unit_amplitude_density`,
  `eval_reaches_the_amplitude_at_the_center`,
  `the_density_peaks_at_the_location_parameter`,
  `the_log_density_is_maximal_at_the_location_parameter`,
  `the_customized_density_is_zero_outside_the_positive_quadrant` - closed-form
  and symmetry invariants of the four densities.

No tier-1 or tier-2 evidence exists for this group: there is no retained C++
output for a `MATH/STATISTICS` fitter and no oracle driver was built.

### Achieved agreement with the published parameters

Measured on the Linux gate node. "Tolerance" is what the C++ class test
asserts; `TEST_REAL_SIMILAR` accepts either an absolute or a relative match, and
OpenMS's defaults are `1e-5` for both.

| Case | Expected | Obtained | Deviation | C++ tolerance |
|---|---|---|---|---|
| `GaussFitter::fit` case 1, `A` | `1.01898275662372` | `1.0189827566259613` | `2.2e-12` rel | default |
| `GaussFitter::fit` case 1, `x0` | `0.300612870901173` | `0.30061287090144295` | `9.0e-13` rel | default |
| `GaussFitter::fit` case 1, `sigma` | `0.136316330927453` | `0.13631633092189843` | `4.1e-11` rel | default |
| `GaussFitter::fit` case 2, `A` | `175011.893006749` | `175011.8930067491` | `6.7e-16` rel | default |
| `GaussFitter::fit` case 2, `x0` | `240.1007246725147` | `240.1007246725147` | exact | default |
| `GaussFitter::fit` case 2, `sigma` | `0.00046642320683761701` | `0.00046642320683761495` | `4.4e-15` rel | default |
| `GaussFitter::eval`, 7 points | see test | all seven | **exact, 7 of 7** | default |
| `GammaDistributionFitter::fit`, `b` | `7.25` | `7.2527011185365495` | `2.7e-3` abs | `0.01` abs |
| `GammaDistributionFitter::fit`, `p` | `3.11` | `3.1168754039854303` | `6.9e-3` abs | `0.01` abs |
| `GumbelDistributionFitter::fit` case 1, `a` | `0.5` | `0.5015861406501456` | `1.6e-3` abs | `0.1` abs |
| `GumbelDistributionFitter::fit` case 1, `b` | `2.0` | `1.9973893576789812` | `2.6e-3` abs | `0.1` abs |
| `GumbelDistributionFitter::fit` case 2, `a` | `1.0` | `0.9955965893616037` | `4.4e-3` abs | `0.1` abs |
| `GumbelDistributionFitter::fit` case 2, `b` | `1.0` | `0.9992689602371769` | `7.3e-4` abs | `0.1` abs |
| `GumbelMaxLikelihoodFitter`, CSV, `a` | `2` | `2.001561295066307` | `1.6e-3` abs | `0.1` abs |
| `GumbelMaxLikelihoodFitter`, CSV, `b` | `0.6` | `0.6230114376387215` | `2.3e-2` abs | `0.1` abs |
| `GumbelMaxLikelihoodFitter`, synthetic, `a` | `2.0` | `1.998676827272269` | `1.3e-3` abs | `0.05` abs |
| `GumbelMaxLikelihoodFitter`, synthetic, `b` | `0.8` | `0.7979922456575415` | `2.0e-3` abs | `0.05` abs |

Nothing fails at the C++ tolerance and nothing was loosened. The Gamma, Gumbel
and maximum-likelihood rows are *not* measurements of the port against the C++:
their expected values are the parameters the data were generated from, rounded
to two or three digits, so the deviations above are dominated by the fit itself.
Only the two `GaussFitter` cases and `GaussFitter::eval` publish the numbers the
C++ actually produced, and those are the rows that read exact to `4e-11`.

The `GaussFitter` rows are the real fidelity measurement. `eval` is now exact
on all seven points, which is the strongest statement available anywhere in this
group: a closed-form expression with no iteration reproduces the C++ bit for
bit. Case 2 starts from a guess already close to its optimum and agrees to
within a few units in the last place. Case 1 starts at `x0 = 3.0` and travels to
`x0 = 0.3` across an order of magnitude more residual evaluations and agrees to
`4e-11`.

**What the `4e-11` is, and what it is not.** Two candidate explanations were
tested rather than asserted:

* It is *not* the `lmpar` association deviation described in §1. Correcting
  that moved case 1 by at most `1.7e-13` relative, more than two orders of
  magnitude short of `4e-11`, and left every other case bit-identical.
* It is *not* the density divisor corrected in §3: `fit` never calls
  `normal_pdf`. Only `eval` does, and `eval` is now exact.

Two candidates remain: the vectorized accumulation order of Eigen's
`squaredNorm` and `dot` (§1), and the `exp` implementation - the Gaussian
residual and its Jacobian call nothing else transcendental. Both produce
last-place differences in the *iterates*, which this case amplifies: it is the
only case that travels a long way across a non-convex surface, and §1's measured
"one ulp in, 810 ulp out" for a deliberately injected one-ulp change is a direct
measurement of that amplification on this exact case. An earlier revision named
`exp` alone as the cause; that was not measured and is not now claimed.

An earlier revision also offered as corroboration that "a Python
re-implementation of the same transcription, run on a different platform and a
different `libm`" reaches the same `f64` values bit for bit. That claim has been
withdrawn. It could not support the conclusion it was attached to: a
re-implementation *of the same transcription* shares every transcription
decision with the Rust, so agreement between them is evidence about language
and nothing about whether the transcription matches Eigen - the very question at
issue. (It would also have had to run on a host whose `libm` differs from the
gate node's to mean what it said, which was never established.) Nothing in this
document now rests on a check that cannot be re-run from this repository.

### Test tolerances

`tests/math_distribution_fitters.rs` asserts, all relative:

| Case | Asserted | Measured deviation | Headroom |
|---|---|---|---|
| `GaussFitter::fit` case 1 | `1e-9` | `4.1e-11` | 24x |
| `GaussFitter::fit` case 2 | `1e-11` | `4.4e-15` | 2,200x |
| `GaussFitter::eval` and `GaussFitResult::eval` | `1e-14` | `0` (bit-exact) | ~45 units in the last place |

Every one of these is tighter than the class test's own `1e-5`, and the
headroom is there for `libm` differences between platforms, not for the port.
An earlier revision of this document said `1e-11` was asserted for `eval`; the
test file asserts `1e-14`, and the doc was simply wrong about its own tests. It
also described the headroom as "roughly two orders of magnitude", which is right
for case 2 and wrong for case 1; the per-row figures above replace it.

`eval` is left at `1e-14` rather than tightened to bit equality: it is bit-exact
on the gate node, but `exp` is not required to be correctly rounded and a
different platform may legitimately move the last bit.

The other fitters use the class tests' own tolerances, because their expected
values are generating parameters rather than C++ output and tightening them
would assert something the source never claimed.

### Checked boundaries

| Boundary | Behaviour |
|---|---|
| More than `MAX_POINTS` observations, or a Jacobian over `MAX_BYTES` | `Error::InvalidValue`, before any allocation |
| Fewer points than parameters | `Error::InvalidValue`; the source reaches the same outcome through `ImproperInputParameters` |
| Non-finite coordinate, sample or weight | `Error::InvalidValue` |
| Negative abscissa for the Gamma fit | `Error::InvalidValue`; `pow(x, p - 1)` is NaN for a negative base |
| Non-finite initial guess, or a zero initial `sigma` / `b` | `Error::InvalidValue`; the first residual would divide by zero |
| `x.len() != w.len()` in `fit_weighted` | `Error::InvalidValue`; the source reads past the end |
| `sigma <= 0` in `eval` or `log_eval_no_normalize` | `Error::InvalidValue`; the source returns NaN from `log` |
| Empty sample set in `fit_weighted` | `Ok`, returning the initial parameters unchanged - the objective is identically zero, so the solver stops at `CosinusTooSmall`. Three class-test sections depend on this |
| Non-finite fitted parameters | `Error::InvalidValue` rather than a silently propagated NaN |
| Every arithmetic path | No `unsafe`, no unchecked indexing, no unchecked integer arithmetic. The internal divisions inside the solver follow IEEE, exactly as Eigen's do; they cannot panic and their denominators are derived from validated input |

---

## 6. Class-test sections

All 25 `START_SECTION` blocks, each with the Rust test and one concrete asserted
value it reproduces. Tests live in `tests/math_distribution_fitters.rs` unless
noted.

### `GaussFitter_test.cpp` (5)

| Section | Rust test | Reproduced value |
|---|---|---|
| `GaussFitter()` | `gauss_fitter_default_construction` | initial guess `(0.06, 3.0, 0.5)` |
| `virtual ~GaussFitter()` | `gauss_fitter_owns_no_resources_to_destroy` | `size_of::<GaussFitter>() == 24`; the C++ section is `NOT_TESTABLE` and only deletes the pointer |
| `GaussFitResult fit(...)` | `gauss_fit_reproduces_both_published_cases` | `A = 1.01898275662372` and `sigma = 0.00046642320683761701` |
| `void setInitialParameters(...)` | `gauss_set_initial_parameters_is_read_back_and_used` | `x0 = 240.1007246725147` from the published guess; `(-1, -1, -1)` read back |
| `static std::vector<double> eval(...)` | `gauss_static_eval_reproduces_the_published_intensities` | `78670.515322697669` at the first m/z |

### `GammaDistributionFitter_test.cpp` (4)

| Section | Rust test | Reproduced value |
|---|---|---|
| `GammaDistributionFitter()` | `gamma_fitter_default_construction` | initial guess `(1.0, 5.0)` |
| `virtual ~GammaDistributionFitter()` | `gamma_fitter_owns_no_resources_to_destroy` | `size_of::<GammaDistributionFitter>() == 16` |
| `GammaDistributionFitResult fit(...)` | `gamma_fit_reproduces_the_published_parameters` | `b = 7.25` at `0.01` absolute |
| `void setInitialParameters(...)` | `gamma_set_initial_parameters_is_read_back_and_used` | `(1.0, 5.0)` read back; `b = 7.25` reached from the `(1.0, 3.0)` guess |

### `GumbelDistributionFitter_test.cpp` (10)

| Section | Rust test | Reproduced value |
|---|---|---|
| `GumbelDistributionFitter()` | `gumbel_fitter_default_construction` | initial guess `(0.25, 0.1)` |
| `virtual ~GumbelDistributionFitter()` | `gumbel_fitter_owns_no_resources_to_destroy` | `size_of::<GumbelDistributionFitter>() == 16` |
| `GumbelDistributionFitResult fit(...)` | `gumbel_fit_reproduces_both_published_cases` | `a = 0.5, b = 2.0` and `a = 1.0, b = 1.0`, both at `0.1` absolute |
| `void setInitialParameters(...)` | `gumbel_set_initial_parameters_is_read_back_and_used` | default-constructed result `(1.0, 2.0)` read back, then `a = 0.5` |
| `GumbelDistributionFitter(const GumbelDistributionFitter&)` | `gumbel_fitter_copy_carries_the_initial_guess` | the copy's guess is `(5.0, 4.0)`; `NOT_TESTABLE` in C++ |
| `GumbelDistributionFitter& operator=(...)` | `gumbel_fitter_assignment_replaces_the_initial_guess` | the target's guess becomes `(3.0, 2.2)`; `NOT_TESTABLE` in C++ |
| `GumbelDistributionFitResult()` | `gumbel_result_default_construction` | `a = 1.0`, `b = 2.0` |
| `GumbelDistributionFitResult(const ...&)` | `gumbel_result_copy_construction` | `a = 5.0`, `b = 4.0` |
| `GumbelDistributionFitResult& operator=(...)` | `gumbel_result_assignment` | `a = 3.0`, `b = 2.2` |
| `MLE` | `gumbel_maximum_likelihood_on_the_published_sample` | 1200 samples read from the fixture; `a = 2`, `b = 0.6` at `0.1` absolute |

### `GumbelMaxLikelihoodFitter_test.cpp` (6)

| Section | Rust test | Reproduced value |
|---|---|---|
| `GumbelMaxLikelihoodFitter()` | `mle_fitter_default_construction` | `fit_weighted(&[], &[])` returns `(0.25, 0.1)` |
| `~GumbelMaxLikelihoodFitter()` | `mle_fitter_owns_no_resources_to_destroy` | `size_of::<GumbelMaxLikelihoodFitter>() == 16` |
| `GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)` | `mle_fitter_construction_from_initial_parameters` | `fit_weighted(&[], &[])` returns `(3.0, 0.5)` |
| `void setInitialParameters(...)` | `mle_set_initial_parameters` | `fit_weighted(&[], &[])` returns `(4.0, 0.7)` |
| `GumbelDistributionFitResult fitWeighted(...)` | `mle_fit_weighted_recovers_the_generating_parameters` | `a = 2.0`, `b = 0.8` at `0.05` absolute, with `b > 0` |
| `[EXTRA] GumbelDistributionFitResult(double, double) + log_eval_no_normalize` | `mle_result_stores_its_parameters_and_evaluates_the_log_density` | `log_eval_no_normalize(2.0) == -1` for `(a, b) = (2.0, 1.0)` |

Unaccounted sections: none.

---

## 7. Source defects found

Recorded here for the integrator; none is worked around silently.

1. **`GumbelDistributionFitter::GumbelDistributionFitResult::eval` is declared
   and never defined** (`GumbelDistributionFitter.h:51`). No definition exists
   anywhere in the SDK, so any C++ caller fails to link. Nothing calls it, which
   is why it has gone unnoticed.
2. **`GumbelDistributionFitter::fitWeighted` is declared and never defined**
   (`GumbelDistributionFitter.h:81`), with a full Doxygen block including an
   `@exception`. Its own class test cannot call it: the `MLE` section of
   `GumbelDistributionFitter_test.cpp` includes `GumbelMaxLikelihoodFitter.h`
   and uses that class instead.
3. **Three headers document a method none of them has, and one translation
   unit calls it.** `GaussFitter.h:29-30`,
   `GumbelDistributionFitter.h:28-29` and `GumbelMaxLikelihoodFitter.h:26-27`
   all say the fitted parameters "can be transformed into a gnuplot formula
   using `getGnuplotFormula()`" (`GaussFitter.h` omits the parentheses). There
   is no such member on any of the three classes; what exists instead is
   formula-building code inside `#ifdef ..._VERBOSE` blocks in the `.cpp` files
   that writes to `std::cout`.

   The stale documentation has a matching stale *caller*, which is worse than
   the comment on its own: `IDDecoyProbability.cpp` calls
   `gdf.getGnuplotFormula()` on a `Math::GammaDistributionFitter` (line 193 and
   line 195) and `gf.getGnuplotFormula()` on a `Math::GaussFitter` (lines 319,
   321 and 328), all inside `#ifdef IDDECOYPROBABILITY_DEBUG`. Defining that
   macro breaks the build of `ANALYSIS/ID`. The only unconditional survivor is
   a commented-out block at `IDDecoyProbability.cpp:311-317`. So the members
   were removed and their callers were disabled rather than updated; the
   getters named `getGnuplotFormula` that still exist in the SDK belong to the
   unrelated `TraceFitter` hierarchy and take four arguments.
4. **`GumbelMaxLikelihoodFitter::fitWeighted` reads past the end of a short
   weight vector** (`GumbelMaxLikelihoodFitter.cpp:58-63`): a second iterator is
   advanced in lockstep with the sample iterator and the two sizes are never
   compared.
5. **The maximum-likelihood formulation is fragile, though correct on the tested
   data.** The objective is placed in `fvec(0)` of a two-element residual vector
   and Levenberg-Marquardt minimizes `||fvec||^2`, so the routine minimizes the
   *square* of the weighted negative log-likelihood. Where that objective is
   strictly positive the stationary points of the square are exactly those of
   the objective and the fit is the maximum-likelihood estimate - verified
   independently for the class-test data by
   `mle_result_is_a_local_minimum_of_the_weighted_log_likelihood`, where the
   objective is `1310.76`. Where the weighted negative log-likelihood can reach
   zero, which a sharply peaked density or small weights make possible, the
   iteration is attracted to that contour instead. The port reproduces the
   formulation and does not correct it.

One further defect is Eigen's, not OpenMS's, and is noted because the port had
to decide what to do about it: `lmpar2` declares `Matrix sdiag(n)` uninitialized
and `qrsolv` reads it to compute `nsing` before the restore that would fill it,
so a `diag[l] == 0` break on the first column reads uninitialized memory. The
driver's `diag` entries are column norms floored at `1.0` and are therefore
always positive, so the break is unreachable from these four fitters; the port
zero-initializes the workspace and says so at the call site.
