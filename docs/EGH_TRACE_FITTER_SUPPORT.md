# EGHTraceFitter

[`src/analysis/feature_finder_picked/egh_trace_fitter.rs`](../src/analysis/feature_finder_picked/egh_trace_fitter.rs)
ports `FEATUREFINDER/EGHTraceFitter.h` and its implementation at core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`: the exponential-Gaussian hybrid
(EGH) retention-time model of Lan and Jorgenson that
`FeatureFinderAlgorithmPicked` fits when `feature:rt_shape` is `asymmetric`.
It is package B5-EGH of the early TOPP bundle
([EARLY_TOPP_WORK_PACKAGES](EARLY_TOPP_WORK_PACKAGES.md), wave 2).

Tests: [`tests/egh_trace_fitter.rs`](../tests/egh_trace_fitter.rs).
Manifest: [`tests/data/egh_trace_fitter_provenance.json`](../tests/data/egh_trace_fitter_provenance.json).
Fixtures: [`tests/data/egh_trace_fitter/`](../tests/data/egh_trace_fitter/).

The header is one Rust file. It implements the `TraceFitter` trait that the
wave-2 scaffold declared in `trace_fitter.rs`, and it uses the shared helpers
package B4-GAUSS added to that file: the start-value estimate
`initial_shape`, the Levenberg-Marquardt driver `optimize`,
`compute_theoretical`, `unable_to_fit` and `stream_number`, and the parameter
defaults. B4-GAUSS owns that file, the Gaussian model and
[TRACE_FITTER_SUPPORT](TRACE_FITTER_SUPPORT.md). The module is not
feature-gated.

**Solver fidelity.** The EGH fit runs through the same `optimize` and
`levenberg_marquardt::minimize` as the Gaussian, and that solver is not yet
bit-faithful to the executed Eigen. The EGH fixtures below agree with the C++
within 1e-9, but that agreement is fixture-specific: on other inputs a fit's
path can depart from Eigen's at its first trial step, and fitted parameters,
status and evaluation count can then differ far beyond 1e-9.
TRACE_FITTER_SUPPORT's "Known gap: solver fidelity beyond the fixtures"
measures this on 79 generated Gaussian inputs and locates the departure inside
`minimize`, not in the functor or the driver configuration, so the caveat
applies to the EGH fits as well. The EGH fits have not been measured beyond
these fixtures. The root cause is in `src/math/fitters/levenberg_marquardt.rs`,
under investigation in lane B3b.

The source marks the class `@experimental`: "Needs further testing on real
data", and its class test exercises the EGH only as a replacement for the
Gaussian, on a symmetric peak. The asymmetric cases below come from a B5 oracle
driver, not from upstream tests.

## API mapping

Every member of the header is listed, including the protected ones and those
inherited from `TraceFitter` that this type implements.

### `EGHTraceFitter::EGHTraceFunctor`

| Source | Rust |
| --- | --- |
| `class EGHTraceFunctor : public TraceFitter::GenericFunctor` | `pub struct EGHTraceFunctor<'a>` |
| `EGHTraceFunctor(int dimensions, const TraceFitter::ModelData* data)` | `EGHTraceFunctor::new(traces: &MassTraces, weighted: bool)`; the dimension is always `EGHTraceFitter::NUM_PARAMS` |
| `~EGHTraceFunctor()` | dropped implicitly |
| `int operator()(const double* x, double* fvec)` | `residuals(&self, x: &[f64], fvec: &mut [f64]) -> Result<()>` |
| `int df(const double* x, double* J)` (column-major) | `jacobian(&self, x: &[f64], jacobian: &mut [f64]) -> Result<()>`, column-major |
| `GenericFunctor::inputs()`, `values()` | `inputs()`, `values()` |
| `const TraceFitter::ModelData* m_data` (protected) | private `traces` and `weighted` fields |

### `EGHTraceFitter`

| Source | Rust |
| --- | --- |
| `EGHTraceFitter()` | `EGHTraceFitter::new()` (`TraceFitterParams::default()`: `max_iteration` 500, unweighted), `Default`, `EGHTraceFitter::with_parameters(TraceFitterParams)` |
| `EGHTraceFitter(const EGHTraceFitter&)`, `operator=` | `Clone` (`clone`, `clone_from`) |
| `~EGHTraceFitter()` | dropped implicitly |
| `void fit(MassTraces&)` | `TraceFitter::fit(&mut self, &MassTraces) -> Result<()>` |
| `double getLowerRTBound() const` | `TraceFitter::lower_rt_bound` |
| `double getUpperRTBound() const` | `TraceFitter::upper_rt_bound` |
| `double getTau() const` | `EGHTraceFitter::tau` (inherent) |
| `double getSigma() const` | `EGHTraceFitter::sigma` (inherent) |
| `double getHeight() const` | `TraceFitter::height` |
| `double getCenter() const` | `TraceFitter::center` |
| `bool checkMaximalRTSpan(double)` | `TraceFitter::check_maximal_rt_span(&self, f64) -> bool` |
| `bool checkMinimalRTSpan(const std::pair<double,double>&, double)` | `TraceFitter::check_minimal_rt_span(&self, (f64, f64), f64) -> bool` |
| `double getValue(double) const` | `TraceFitter::value` |
| `double getArea()` | `TraceFitter::area(&self)` |
| `double getFWHM() const` | `TraceFitter::fwhm` |
| `std::string getGnuplotFormula(const MassTrace&, char, double, double)` | `TraceFitter::gnuplot_formula(&self, &MassTrace, char, f64, f64) -> String` |
| protected `double apex_rt_, height_, sigma_, tau_` | private fields, read through `center`, `height`, `sigma`, `tau` |
| protected `std::pair<double,double> sigma_5_bound_` | private field, read through `lower_rt_bound` and `upper_rt_bound` |
| protected `double region_rt_span_` | private field, used by `check_maximal_rt_span` and exposed on the start point as `EGHInitialParameters::region_rt_span` |
| protected `static const double EPSILON_COEFS_[]` | `EGHTraceFitter::EPSILON_COEFS: [f64; 7]` (public) |
| protected `static const Size NUM_PARAMS_` | `EGHTraceFitter::NUM_PARAMS: usize` (public) |
| protected `getAlphaBoundaries_(double alpha) const` | `EGHTraceFitter::alpha_boundaries(&self, f64) -> (f64, f64)` (public) |
| protected `getOptimizedParameters_(const std::vector<double>&)` | `EGHTraceFitter::set_optimized_parameters(&mut self, [f64; 4])` (public) |
| protected `setInitialParameters_(MassTraces&)` | `EGHTraceFitter::initial_parameters(&MassTraces) -> Result<EGHInitialParameters>` (public, associated), on the shared `trace_fitter::initial_shape` with `ProfileSmoothing::Always` |
| protected `updateMembers_()` | `TraceFitter::set_parameters`; the source only calls the base implementation |
| `@htmlinclude OpenMS_EGHTraceFitter.parameters` | the two `TraceFitter` parameters, `TraceFitterParams`; their `Param` mapping is package B4-GAUSS's |

### Inherited from `TraceFitter`

| Source | Rust |
| --- | --- |
| `computeTheoretical(const MassTrace&, Size) const` | `TraceFitter::compute_theoretical(&self, &MassTrace, usize) -> Result<f64>`, delegating to the shared `trace_fitter::compute_theoretical` |
| `getParameters`, `setParameters` (`DefaultParamHandler`) | `TraceFitter::parameters`, `TraceFitter::set_parameters` |
| protected `optimize_(std::vector<double>&, GenericFunctor&)` | the shared `trace_fitter::optimize` |
| protected `SignedSize max_iterations_`, `bool weighted_` | `TraceFitterParams::max_iteration: i64`, `weighted: bool`; defaults and `Param` mapping in `trace_fitter.rs` |
| protected `struct ModelData` | internal to `EGHTraceFunctor` |

The three protected hooks are public here: `initial_parameters` lets a caller
inspect a fit's start point, and `set_optimized_parameters` rebuilds a model
from stored parameters, for example a feature's `EGH_height`, `EGH_tau` and
`EGH_sigma` meta values with its retention time. `EGHInitialParameters` is a
Rust record for the members `setInitialParameters_` writes.

`FeatureFinderAlgorithmPicked` reads `getTau`, `getHeight` and `getSigma`
through `std::dynamic_pointer_cast<EGHTraceFitter>` on the fitter
`chooseTraceFitter_` returns. The trait has no downcast, so the later algorithm
package needs to hold the concrete type (for example in an enum of the two
fitters) where it writes those meta values.

## Preserved source conventions

- **Operand order.** Every expression keeps the source's order and grouping:
  `2 * sigma * sigma` is `(2 sigma) sigma`, `-1 / log(alpha) * (B - A)` divides
  first, `-0.5 / log(alpha) * B * A` multiplies left to right, and the Jacobian
  products are formed left to right before the division by `d * d`. The
  source's `-1 * (L * tau_)` is written as a negation, which gives the same bits.
- **Residual.** `(baseline + theo * H * exp(-t^2 / d) - I) * w` for `d > 0`,
  `(0 - I) * w` otherwise: the baseline is dropped where the denominator is not
  positive. The intensity is the `f32` peak intensity promoted to `f64`.
- **Signed and absolute sigma.** The residual uses the signed `sigma`; the
  Jacobian uses `|sigma|` ("must be non-negative!"), so for a negative `sigma`
  its sigma column has the opposite sign of the true derivative.
- **Weighting.** `w = theoretical_int` when `weighted`, `1` otherwise, applied
  to residuals and Jacobian rows alike.
- **Row order.** Traces in order, then each trace's peaks in order.
- **Start point.** The shared `initial_shape` with `ProfileSmoothing::Always`:
  the intensity profile of `MassTraces::intensity_profile`, a zero-padded
  running sum over five entries divided by `5`, seeded as
  `std::accumulate` seeds it (`0.0 + totals[2] + totals[3]`); the first strict
  maximum; the half-height walks with strict `>`; `alpha = (left + right) * 0.5
  / height`. There is no guard for short profiles and none for `alpha >= 1`,
  unlike the Gaussian model: `alpha = 1` gives `sigma = NaN` and `tau = -inf`
  or NaN, `alpha > 1` a NaN `sigma`. A `tau` of exactly zero (also `-0.0`)
  becomes `f64::EPSILON`.
- **Status handling.** Through the shared `optimize`, every
  Levenberg-Marquardt status after `ImproperInputParameters` is accepted with
  the solver's vector, including
  `TooManyFunctionEvaluation` and `CosinusTooSmall` at the start point, which is
  how a NaN start ends: its Jacobian is all zero.
- **Error wording.** `UnableToFit-FinalSet` with "Skipping feature, we always
  expect N>=p" for fewer than four peaks, and "Could not fit the gaussian to the
  data: Error 0" for a non-positive `max_iteration`, the source's text even for
  the EGH.
- **Bounds.** `sigma_5_bound_` at relative height `0.043937`, the FWHM from
  relative height `0.5`, `std::min`/`std::max` semantics for the ordered pair:
  the first operand wins unless the second is strictly smaller (larger), so a
  NaN second operand yields the first and `+0.0` is kept against `-0.0`.
- **Area.** `phi = atan(|tau| / |sigma|)`, `epsilon` accumulated from
  `EPSILON_COEFS[0]` in ascending powers of `phi`, `H * (|sigma| * 0.6266571 +
  |tau|) * epsilon`. `sigma = tau = 0` gives NaN.
- **Span checks.** `check_minimal_rt_span` is `rt span < min_rt_span * bound
  span`, `check_maximal_rt_span` is `bound span > max_rt_span * region span`;
  `true` reports the violation that `checkFeatureQuality_` rejects, as the trait
  documents.
- **Gnuplot formula.** The source's exact text, with every number written as a
  default C++ stream writes a `double` (precision 6, `%g` style, `-0`, `nan`,
  `inf`), through the shared `trace_fitter::stream_number`.
- **Serial.** The source is serial, and so is the port.

## Native differences

1. **Atomic failure.** A failed `fit` leaves the fitter unchanged. The source has
   already written the start point into `height_`, `apex_rt_`, `sigma_`, `tau_`
   and `region_rt_span_` when `optimize_` throws, and keeps the previous
   `sigma_5_bound_`. `FeatureFinderAlgorithmPicked` does not catch the exception
   inside its seed loop, so no caller observes that partial state.
   `GaussTraceFitter` keeps its start values after a refused fit, as its
   source does, so the two `TraceFitter` implementations differ here; the
   trait leaves the post-error state to each implementation, and one policy
   for both is the integrator's decision.
2. **Traces without peaks.** `fit` returns `UnableToFit-FinalSet: Skipping
   feature, we always expect N>=p` and `initial_parameters` returns
   `initial_shape`'s `Error::InvalidValue`. The source's
   `setInitialParameters_` reads `smoothed[0]` of an empty vector, which is
   undefined behaviour, before `optimize_` would throw.
3. **NaN retention times.** A NaN retention time that meets a profile entry
   while the profile is merged gives `MassTraces::intensity_profile`'s
   `Error::InvalidValue`; the source loop never terminates.
4. **Checked indices and buffers.** `compute_theoretical` refuses an
   out-of-range peak index, and `EGHTraceFunctor::residuals` and `jacobian`
   refuse mis-sized slices. The source indexes and writes through raw pointers
   without a check.
5. **Initial values.** The model values are `0.0` before the first fit;
   `set_optimized_parameters` on a new fitter leaves the region span `0.0`. The
   source leaves both uninitialised.
6. **Protected members exposed.** `initial_parameters`,
   `set_optimized_parameters`, `alpha_boundaries`, `EPSILON_COEFS` and
   `NUM_PARAMS` are public.
7. **No debug log.** `setInitialParameters_` logs its intermediate values at
   debug level; the port logs nothing.
8. **Mathematical functions.** `exp`, `log`, `sqrt` and `atan` come from the
   `libm` crate, so results are the same on every platform up to the sign and
   payload of a NaN. The libm crate's `exp` is not correctly rounded, and
   single evaluations differ from the oracle's Apple libm by one or two units
   in the last place on about 6% of the functor values; Levenberg-Marquardt
   results differ by at most 2.0e-11 relative. See the evidence below.
   `GaussTraceFitter` calls the platform `exp` and `log` instead
   (TRACE_FITTER_SUPPORT, native difference 1 and "`exp` and `log` across
   platforms"), which asks that both fitters make the same choice. That choice
   is the integrator's; this module keeps the `libm` crate until it is made.
9. **Function name.** The gnuplot function name is a Rust `char`; a non-ASCII
   character is written as its UTF-8 bytes, where the source writes one byte.
10. **Work ceilings of the shared driver.** `fit` inherits `optimize`'s
    ceilings (TRACE_FITTER_SUPPORT, native difference 5): the solver's point
    and byte ceilings are checked before anything is evaluated, and the budget
    passed to the solver is at most `MAX_RESIDUAL_WORK / values` (2^30 residual
    evaluations in total). A fit that exhausts that ceiling before its
    configured `max_iteration` is refused with `Error::InvalidValue`, and the
    fitter is unchanged, where the source would continue or accept. Every EGH
    fixture ends far below the ceiling, so its results are those of the
    uncapped budget.

## Checked boundaries and evidence

### Evidence tiers

- **Tier 1, executed differential.** Two product-SDK oracles, both with inputs
  recorded bit for bit (2799 rows over 28 cases):
  - `egh_trace_fitter_c2.tsv`: the EGH library records of the C2 class-level
    oracle (`../oracle/featurefinder-picked`), extracted by
    `../oracle/egh-trace-fitter/extract_c2.py`. Cases: the class-test traces at
    theoretical intensities 0.8/0.2 and 0.4/0.6, each weighted and unweighted,
    and the degenerate inputs `flat3`, `short3` and `values_lt_inputs`.
  - `egh_trace_fitter_oracle.tsv`: `../oracle/egh-trace-fitter/driver.cpp`.
    Cases: `tailing_3traces` and its weighted twin (three traces on
    staggered retention-time ranges, baseline 11.25, generated with tau 2.4),
    `fronting_2traces_interleaved` (tau -1.8, grids half a step apart),
    `single_trace_4peaks` (values equal to parameters), `single_trace_3peaks`
    (refused), `alpha_above_one_baseline`, `apex_at_last_scan` (budget
    exhausted), `max_iteration_0`, `max_iteration_minus_1`, the defaults, and
    eleven parameter vectors (tailing, negative sigma and tau, tau 0, all
    widths 0, sigma 0, NaN height, large fronting, `tau/sigma = 1e6`, negative
    height, signed zeros, infinite tau) with bounds at seven alphas including
    `0`, `1` and negative.
- **Tier 3, source review.** The 18 `START_SECTION`s of `EGHTraceFitter_test.cpp`.
- **Tier 4, Rust-only.** Atomic failure, traces without peaks, the NaN merge
  error, index and buffer checks, the trait object, `Clone` and `Default`.

### Class test

| Section | Rust test | Expectation |
| --- | --- | --- |
| `EGHTraceFitter()` | `section_default_constructor` | constructed; defaults 500, unweighted |
| `~EGHTraceFitter()` | `section_destructor` | dropped |
| copy constructor | `section_copy_constructor` | center, height, bounds equal |
| `operator=` | `section_assignment_operator` | center, height, bounds equal |
| `fit` | `section_fit` | x0 680.1, H 10; weighted H 10; weighted 0.4/0.6 H 6.0825 |
| `getLowerRTBound` | `section_get_lower_rt_bound` | 680.1 - 3.75 |
| `getUpperRTBound` | `section_get_upper_rt_bound` | 680.1 + 3.75 |
| `getHeight` | `section_get_height` | 10 |
| `getCenter` | `section_get_center` | 680.1 |
| `getTau` | `section_get_tau` | 0 (absolute 1e-5) |
| `getSigma` | `section_get_sigma` | 1.5 |
| `getValue` | `section_get_value` | 10 at 680.1 |
| `computeTheoretical` | `section_compute_theoretical` | 8 |
| `checkMaximalRTSpan` | `section_check_maximal_rt_span` | false, then true after `- 0.1` |
| `checkMinimalRTSpan` | `section_check_minimal_rt_span` | false at 0.5, true at 1.0 |
| `getArea` | `section_get_area` | `sqrt(2 pi) * 1.5 * 10` |
| `getGnuplotFormula` | `section_get_gnuplot_formula` | prefix, three substrings, suffix |
| `getFWHM` | `section_get_fwhm` | 3.53223007592464 |

`TEST_REAL_SIMILAR` is `isRealSimilar` with absolute 1e-5 and ratio 1 + 1e-5.
The class test's first `checkMaximalRTSpan` argument makes the limit exactly
`5 * 1.5 = 7.5`. Relative height `0.043937` lies `2.4999994` sigma from the
apex, so the fitted bound span is `7.4999982` and `false` holds with a margin of
`1.8e-6`, for the C++ and for the port.

`c2_acceptance_values` checks the package acceptance values in decimal within
1e-9 relative: H 9.9999999742433339, sigma 1.5000000035584744, FWHM
3.5322300759260088, area 37.599425992354142, weighted 0.4/0.6 H
6.0824742044607856.

### Comparison classes

| Class | What | Tolerance |
| --- | --- | --- |
| Exact | booleans, integers, `UnableToFit` messages, gnuplot formulas of a model set to the oracle's parameters, getters after `set_optimized_parameters`, the budget-boundary pattern | equality |
| Single evaluation | functor residuals and Jacobians at recorded vectors, start points, bounds, FWHM, area, values, `computeTheoretical` and alpha boundaries of a model set to the oracle's parameters | `|a - e| <= 1e-14 * max(|a|, |e|)`; a residual on `max(|a|, |e|, |I w|)` |
| Levenberg-Marquardt | fitted parameters, their bounds, FWHM, area and `computeTheoretical`, and every recorded budget of the sweep | `|a - e| <= 1e-9 * max(|a|, |e|)`; a `tau` at rounding noise on `|sigma|` |

NaN matches NaN whatever its sign bit. The queries of the fitted model are
evaluated on a clone of the Rust fit set to the oracle's parameters, so they
keep the Rust fit's region span, a plain difference of recorded retention
times, and compare a single evaluation rather than accumulated
Levenberg-Marquardt differences.

`tau` of a symmetric peak stays at rounding noise (about `4e-15` on the class
test), whose relative value carries no information; the model sees `tau` only
as `tau * t` next to `2 sigma^2`. So where both values of `tau` lie below
`1e-9 * max(|sigma|)`, `tau` is compared on the scale of `sigma`; everywhere
else it gets the plain relative bound. The sigma scale applies to 24 of the 594
`tau` comparisons, all with `|tau| <= 8.5e-15` and differences of at most
`5.6e-16` absolute.

### Measurements

On Linux x86-64 (IBMI node dax), stable and 1.85.0:

| Test | Comparisons | Bit-identical or equal |
| --- | --- | --- |
| `oracle_initial_parameters` | 80 | 80 |
| `oracle_functor_residuals_and_jacobians` | 11023 | 10380 |
| `oracle_parameter_vectors` | 515 | 510 |
| `oracle_fit_and_queries` | 1260 | 1003 |
| `oracle_budget_boundaries` | 6828 | 4680 |

- The differing single evaluations are one or two units in the last place of an
  `exp` term. On the 3600 distinct `exp` inputs of the functor rows, Apple's
  `exp` in the oracle returned the correctly rounded value 3594 times (checked
  against Python `decimal` at 80 digits), while the libm crate's `exp` differs
  on about 6% of the functor values. A temporary run with glibc's `exp`
  through `f64::exp` left 14 functor values different, but the fits still
  differed in their last places, so the Levenberg-Marquardt transcription
  itself is not bit-identical to Eigen on these problems. Bit identity with the
  C++ therefore needs both a correctly rounded `exp` and a bit-identical solver,
  and is not claimed.
- **Platform-matched check** (review of `45aa00c`, repeated on the rebased
  module). On macOS arm64, the oracle platform, a scratch copy replaced the
  libm crate's `exp`, `log` and `atan` with `f64::exp`, `f64::ln` and
  `f64::atan` (Apple libm). Start points then matched 80 of 80 bit for bit,
  functor values 11023 of 11023 and parameter vectors 515 of 515, so the
  residual, Jacobian, start-point and query transcriptions are exact. Fits
  stayed at 1125 of 1260 and the budget sweep at 4755 of 6828, which confirms
  on the oracle platform itself that the Levenberg-Marquardt path is not
  bit-identical to Eigen here (see "Solver fidelity" above).
- **Acceptance criterion 4** of the bundle plan ("residual and Jacobian equal
  C2 at recorded vectors") is therefore met within 1e-14 relative with the
  mandated `libm` crate (10380 of 11023 bitwise), and bitwise only with the
  platform library on the oracle platform. Whether bitwise equality is
  required is the lead decision the package asked for; it is not yet
  recorded.
- Largest Levenberg-Marquardt difference: 2.0e-11 relative, `tau` on the
  `apex_at_last_scan` budget sweep at budget 240; the case never converges and
  spends every budget while `sigma` and `tau` grow. After all 500 evaluations
  its parameters differ by at most 1.4e-11 (`tau`) and its derived bounds,
  FWHM and area by at most 1.5e-11. Elsewhere the parameters differ by at most
  3.9e-16 relative, `tau` at rounding noise by up to 5.6e-16 absolute, the
  derived bounds, FWHM and area by at most 4.1e-16, and `computeTheoretical` of
  the Rust fit far in the tails (values near `1e-27`) by up to 1.2e-12
  relative.
- The Rust results are the same on every platform up to the sign and payload
  of a NaN: the `libm` crate (`=0.2.16`) is pure Rust, Rust does not contract
  floating-point operations, and the solver uses only IEEE basic operations and
  `sqrt`. The review hashed every compared Rust value on macOS arm64 and Linux
  x86-64: functor, fit and budget values agree, while start-point and
  parameter-vector hashes differ only in the sign bit of NaN results (x86-64
  produces negative default NaNs). After the rebase on B4 every value the test
  file computes was dumped on macOS arm64 and on dax: the two dumps agree
  except for the sign bit of NaN results. The tests treat every NaN as equal,
  and the gnuplot formula writes any NaN as `nan`. A change of the `libm` version or of
  the Levenberg-Marquardt backend (packages B3 and B3b) needs this test re-run
  and, if it moves, a new measurement, not a wider tolerance.

### Evaluation budget

`max_iteration` is Eigen's `maxfev`. For every swept case and every budget
`n = 1 .. 500`, the Rust fit at `n` equals the Rust fit at 500 bit for bit
exactly where the C++ fit at `n` equals the C++ fit at 500:

| Case | Identical from |
| --- | --- |
| class test, 0.8/0.2 and 0.4/0.6, weighted and unweighted | 4 |
| `tailing_3traces` | 6 |
| `tailing_3traces_weighted` | 7 |
| `fronting_2traces_interleaved` | 7 |
| `single_trace_4peaks` | 38 |
| `apex_at_last_scan` | 500 (every smaller budget differs) |

The recorded budgets below and one past each boundary (all 500 for
`apex_at_last_scan`) also match the C++ parameters within 1e-9. This pins the
evaluation accounting and termination order of the Levenberg-Marquardt backend
on the EGH problems; package B3's crate adapter must keep it.

### Degenerate inputs

| Input | C++ and port |
| --- | --- |
| `flat3` (alpha 1) | start `sigma` NaN, `tau` -inf; the fit is accepted unchanged |
| `short3` (three scans, alpha 1) | start `sigma` NaN, `tau` -inf; accepted unchanged |
| `alpha_above_one_baseline` | alpha 1.2: `sigma` NaN, `tau` -2.74241; accepted, NaN bounds, FWHM and area, all values 0 |
| `values_lt_inputs` (2 peaks), `single_trace_3peaks` | `UnableToFit-FinalSet: Skipping feature, we always expect N>=p` |
| `max_iteration` 0 and -1 | `UnableToFit-FinalSet: Could not fit the gaussian to the data: Error 0` |
| traces without peaks | C++ undefined behaviour (not run); port `UnableToFit-FinalSet` |

### Resources and performance

A fit allocates the intensity profile and the smoothed profile once, each
proportional to the number of distinct retention times (`initial_shape` reads
the zero padding through an index function), plus `optimize`'s copy of the
parameter vector and the solver's buffers; the residual and Jacobian loops
allocate nothing and walk the traces once per evaluation. The profile's peak
and merge-step ceilings and the solver's point and byte ceilings are checked
before those allocations, and `optimize` caps the total residual work.

## C++ issue candidates

1. **No `alpha >= 1` guard in `setInitialParameters_`** (executed). A flat or
   single-scan profile, or a baseline above the smoothed half-height edges,
   gives a NaN `sigma` and a meaningless `tau`, and `fit` succeeds with those
   values (`CosinusTooSmall` at the start point), with NaN bounds, FWHM and
   area. `GaussTraceFitter` guards the same case. Port: reproduced.
2. **Empty traces** (source review, undefined behaviour). `setInitialParameters_`
   reads `smoothed[0]` of an empty vector. Port: `UnableToFit-FinalSet`.
3. **Baseline and sigma sign in the functor** (source review and executed rows).
   Where the denominator is not positive the residual model is `0`, not the
   baseline, while `getGnuplotFormula` adds the baseline; `df` uses `|sigma|`
   and the residual the signed `sigma`, so the sigma column has the wrong sign
   for `sigma < 0`. Port: reproduced.
4. **Wording** (executed). The `UnableToFit` text says "Could not fit the
   gaussian to the data" for the EGH model. Port: reproduced.

## Ledger notes

- Suggested status for `FEATUREFINDER/EGHTraceFitter.h`: `complete`. The
  module is rebased on B4's `trace_fitter::optimize`, defaults and helpers, and
  every value `tests/egh_trace_fitter.rs` computes is bit-identical to the
  pre-rebase module (the only change is the error text of `initial_parameters`
  on traces without peaks, now `initial_shape`'s). Record the same known gap
  as for `GaussTraceFitter`: the fit is not bit-faithful to the executed
  library beyond the fixtures (solver fidelity, lane B3b).
- No module edge of its own: the module reaches `crate::math` only through
  `trace_fitter::optimize` (B4's `analysis -> math`), and `crate::format` only
  through `trace_fitter::stream_number`.
