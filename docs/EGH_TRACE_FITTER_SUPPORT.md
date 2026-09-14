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
wave-2 scaffold declared in `trace_fitter.rs` (package B4-GAUSS owns that file
and the Gaussian model). The module is not feature-gated.

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
| `EGHTraceFitter()` | `EGHTraceFitter::new()` (defaults `max_iteration` 500, unweighted), `Default`, `EGHTraceFitter::with_parameters(TraceFitterParams)` |
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
| protected `setInitialParameters_(MassTraces&)` | `EGHTraceFitter::initial_parameters(&MassTraces) -> Result<EGHInitialParameters>` (public, associated) |
| protected `updateMembers_()` | `TraceFitter::set_parameters`; the source only calls the base implementation |
| `@htmlinclude OpenMS_EGHTraceFitter.parameters` | the two `TraceFitter` parameters, `TraceFitterParams`; their `Param` mapping is package B4-GAUSS's |

### Inherited from `TraceFitter`

| Source | Rust |
| --- | --- |
| `computeTheoretical(const MassTrace&, Size) const` | `TraceFitter::compute_theoretical(&self, &MassTrace, usize) -> Result<f64>` |
| `getParameters`, `setParameters` (`DefaultParamHandler`) | `TraceFitter::parameters`, `TraceFitter::set_parameters` |
| protected `optimize_(std::vector<double>&, GenericFunctor&)` | private `optimize_stand_in`, to be replaced by B4's `trace_fitter::optimize` |
| protected `SignedSize max_iterations_`, `bool weighted_` | `TraceFitterParams::max_iteration: i64`, `weighted: bool` |
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
- **Start point.** The intensity profile of `MassTraces::intensity_profile`,
  a zero-padded running sum over five entries divided by `5`, seeded as
  `std::accumulate` seeds it (`0.0 + totals[2] + totals[3]`); the first strict
  maximum; the half-height walks with strict `>`; `alpha = (left + right) * 0.5
  / height`. There is no guard for short profiles and none for `alpha >= 1`,
  unlike the Gaussian model: `alpha = 1` gives `sigma = NaN` and `tau = -inf`
  or NaN, `alpha > 1` a NaN `sigma`. A `tau` of exactly zero (also `-0.0`)
  becomes `f64::EPSILON`.
- **Status handling.** Every Levenberg-Marquardt status after
  `ImproperInputParameters` is accepted with the solver's vector, including
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
  `inf`), through `format::file_info::text_format::ostream_g`.
- **Serial.** The source is serial, and so is the port.

## Native differences

1. **Atomic failure.** A failed `fit` leaves the fitter unchanged. The source has
   already written the start point into `height_`, `apex_rt_`, `sigma_`, `tau_`
   and `region_rt_span_` when `optimize_` throws, and keeps the previous
   `sigma_5_bound_`. `FeatureFinderAlgorithmPicked` does not catch the exception
   inside its seed loop, so no caller observes that partial state.
2. **Traces without peaks.** `fit` returns `UnableToFit-FinalSet: Skipping
   feature, we always expect N>=p` and `initial_parameters` returns
   `Error::InvalidValue`. The source's `setInitialParameters_` reads
   `smoothed[0]` of an empty vector, which is undefined behaviour, before
   `optimize_` would throw.
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
   `libm` crate, so results are the same on every platform. The libm crate's
   `exp` is not correctly rounded, and single evaluations differ from the
   oracle's Apple libm by one or two units in the last place on about 6% of the
   functor values; Levenberg-Marquardt results differ by at most 2e-11
   relative. See the evidence below.
9. **Function name.** The gnuplot function name is a Rust `char`; a non-ASCII
   character is written as its UTF-8 bytes, where the source writes one byte.
10. **Temporary driver (until B4 merges).** `fit` calls the crate's
    Levenberg-Marquardt `minimize` through the private `optimize_stand_in`,
    which has the signature fixed in `trace_fitter.rs` for B4's `optimize` and
    the semantics of `TraceFitter::optimize_`: refuse fewer residuals than
    parameters, preflight the Jacobian size, `max_fev = max_iteration` (zero or
    below becomes a zero budget, which the solver rejects as improper input),
    refuse a status of `ImproperInputParameters` or below. The TraceFitter
    defaults are a private copy. Both go when this module is rebased on B4.

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
| Levenberg-Marquardt | fitted parameters, their bounds, FWHM, area and `computeTheoretical`, and every recorded budget of the sweep | `|a - e| <= 1e-9 * max(|a|, |e|)`; `tau` on `max(|a|, |e|, |sigma|)` |

NaN matches NaN whatever its sign bit. The queries of the fitted model are
evaluated on a clone of the Rust fit set to the oracle's parameters, so they
keep the Rust fit's region span, a plain difference of recorded retention
times, and compare a single evaluation rather than accumulated
Levenberg-Marquardt differences.

`tau` of a symmetric peak stays at rounding noise (about `4e-15` on the class
test), whose relative value carries no information; the model sees `tau` only
as `tau * t` next to `2 sigma^2`, so it is compared on the scale of `sigma`.

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
- Largest Levenberg-Marquardt difference: 2.0e-11 relative, on
  `apex_at_last_scan`, which spends all 500 evaluations while `sigma` and `tau`
  grow. Elsewhere the parameters differ by at most 3.9e-16 relative, `tau` at
  rounding noise by up to 5.6e-16 absolute, the derived bounds, FWHM and area by
  at most 7.7e-16, and `computeTheoretical` of the Rust fit far in the tails
  (values near `1e-27`) by up to 1.2e-12 relative.
- The Rust results themselves do not depend on the platform: the `libm` crate
  (`=0.2.16`) is pure Rust, Rust does not contract floating-point operations,
  and the solver uses only IEEE basic operations and `sqrt`. A change of the
  `libm` version or of the Levenberg-Marquardt backend (package B3) needs this
  test re-run and, if it moves, a new measurement, not a wider tolerance.

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

A fit allocates the intensity profile, the padded totals and the smoothed
profile once, each proportional to the number of distinct retention times, and
the solver's buffers; the residual and Jacobian loops allocate nothing and walk
the traces once per evaluation. The profile's peak and merge-step ceilings and
the solver's point and byte ceilings are checked before those allocations.

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

- Suggested status for `FEATUREFINDER/EGHTraceFitter.h`: `partial` while the
  private Levenberg-Marquardt stand-in and default copy remain; `complete` once
  the module is rebased on B4's `trace_fitter::optimize` and defaults and this
  test passes unchanged there.
- New module edge: `analysis -> math` (the stand-in's `minimize`), which the
  scaffold lists for B4 and which closes no cycle. `analysis -> format` (for
  `ostream_g`) already exists.
