# Trace fitters: TraceFitter and GaussTraceFitter

[`src/analysis/feature_finder_picked/trace_fitter.rs`](../src/analysis/feature_finder_picked/trace_fitter.rs)
ports `FEATUREFINDER/TraceFitter.h` and `TraceFitter.cpp`, and
[`src/analysis/feature_finder_picked/gauss_trace_fitter.rs`](../src/analysis/feature_finder_picked/gauss_trace_fitter.rs)
ports `FEATUREFINDER/GaussTraceFitter.h` and `GaussTraceFitter.cpp`, at core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. `FeatureFinderAlgorithmPicked` fits
the Gaussian to the mass traces of every feature candidate when
`feature:rt_shape` is `symmetric` (the default and FeatureFinderCentroided_1's
setting). This is package B4-GAUSS of the early TOPP bundle
([`EARLY_TOPP_WORK_PACKAGES.md`](EARLY_TOPP_WORK_PACKAGES.md)); the EGH model
(`EGHTraceFitter.h`) is package B5-EGH and builds on the same trait.

Tests: [`tests/trace_fitter.rs`](../tests/trace_fitter.rs),
[`tests/gauss_trace_fitter.rs`](../tests/gauss_trace_fitter.rs), and the unit
tests of the work ceiling in `trace_fitter.rs`.
Fixtures: [`tests/data/gauss_trace_fitter/`](../tests/data/gauss_trace_fitter/).
Manifest: [`tests/data/gauss_trace_fitter_provenance.json`](../tests/data/gauss_trace_fitter_provenance.json).

Nothing is feature-gated. The modules name `crate::math` (the
Levenberg-Marquardt solver) and `crate::param` (the defaults), two new acyclic
`analysis` edges, plus `crate::format` (the `std::ostream` number text) and
`crate::analysis`, which `analysis` already had. The fitters are serial, as the
source is.

## API mapping

### `TraceFitter`

The abstract class becomes the object-safe trait `TraceFitter` and the record
`TraceFitterParams`. The trait's signatures and the record's fields are the
wave-2 integration contract declared by the integrator; this package added the
implementations and helpers below them.

| Source | Rust |
| --- | --- |
| `class TraceFitter : public DefaultParamHandler` | `trait TraceFitter` |
| `TraceFitter()` (registers `max_iteration` 500 and `weighted` `"false"`, both `advanced`) | `TraceFitterParams::default()`; `TraceFitterParams::defaults() -> Result<Param>` for the `Param` with descriptions, tags and valid strings |
| `TraceFitter(const TraceFitter&)`, `operator=` | `Clone` on each implementor |
| `~TraceFitter()` | `Drop` (implicit) |
| `getParameters` / `setParameters` + `updateMembers_` (inherited) | `TraceFitter::parameters` / `TraceFitter::set_parameters(TraceFitterParams)`; `TraceFitterParams::from_param(&Param) -> Result<(Self, Vec<String>)>` and `TraceFitterParams::to_param` convert |
| `getDefaults` (inherited) | `TraceFitterParams::defaults` |
| `virtual void fit(MassTraces&) = 0` | `fn fit(&mut self, &MassTraces) -> Result<()>` |
| `getLowerRTBound`, `getUpperRTBound`, `getHeight`, `getCenter`, `getFWHM` | `lower_rt_bound`, `upper_rt_bound`, `height`, `center`, `fwhm` |
| `virtual double getValue(double rt) const = 0` | `fn value(&self, rt: f64) -> f64` |
| `double computeTheoretical(const MassTrace&, Size k) const` (non-virtual) | required `fn compute_theoretical(&self, &MassTrace, usize) -> Result<f64>`; shared body `trace_fitter::compute_theoretical` |
| `virtual bool checkMinimalRTSpan(const std::pair<double,double>&, double) = 0` | `fn check_minimal_rt_span(&self, (f64, f64), f64) -> bool` |
| `virtual bool checkMaximalRTSpan(double) = 0` | `fn check_maximal_rt_span(&self, f64) -> bool` |
| `virtual double getArea() = 0` | `fn area(&self) -> f64` |
| `virtual std::string getGnuplotFormula(const MassTrace&, char, double, double) = 0` | `fn gnuplot_formula(&self, &MassTrace, char, f64, f64) -> String`; numbers through `trace_fitter::stream_number` |
| `class GenericFunctor` (`inputs`, `values`, `operator()`, `df`, `m_inputs`, `m_values`) | the closures passed to `optimize`; `inputs()` is `x.len()`, `values()` the `values` argument |
| `GenericFunctorEigenAdapter` (file-local) | the `jacobian` wrapper inside `optimize_with_status`, which returns 0 consumed evaluations |
| `struct ModelData { traces_ptr, weighted }` (protected) | fields of each fitter's functor |
| `virtual void getOptimizedParameters_(const std::vector<double>&) = 0` (protected) | an inherent method of each fitter (`GaussTraceFitter::set_optimized_parameters`) |
| `void optimize_(std::vector<double>&, GenericFunctor&)` (protected) | `optimize(x, values, residual, jacobian, &TraceFitterParams) -> Result<()>`; `optimize_with_status` returns the accepted `LmStatus` |
| `Exception::UnableToFit` named `UnableToFit-FinalSet` | `Error::InvalidValue("UnableToFit-FinalSet: <message>")` from `unable_to_fit`; constants `UNABLE_TO_FIT_FINAL_SET`, `FEWER_RESIDUALS_THAN_PARAMETERS` |
| `SignedSize max_iterations_` (protected) | `TraceFitterParams::max_iteration: i64` |
| `bool weighted_` (protected) | `TraceFitterParams::weighted: bool` |
| (native) | `TraceFitterParams::DEFAULT_MAX_ITERATION`, `HANDLER_NAME`, `MAX_ITERATION_DESCRIPTION`, `WEIGHTED_DESCRIPTION`; `MAX_RESIDUAL_WORK`; `initial_shape`, `InitialShape`, `ProfileSmoothing` (the start-value steps both subclasses share) |

### `GaussTraceFitter`

| Source | Rust |
| --- | --- |
| `GaussTraceFitter()` | `GaussTraceFitter::new`, `Default`; `with_parameters(TraceFitterParams)` |
| `GaussTraceFitter(const GaussTraceFitter&)`, `operator=` | `Clone` (see native differences) |
| `~GaussTraceFitter()` | `Drop` (implicit) |
| `fit`, `getLowerRTBound`, `getUpperRTBound`, `getHeight`, `getCenter`, `getFWHM`, `checkMaximalRTSpan`, `checkMinimalRTSpan`, `getValue`, `getArea`, `getGnuplotFormula` | the `TraceFitter` impl |
| `double getSigma() const` | `GaussTraceFitter::sigma` |
| `double sigma_`, `x0_`, `height_` (protected) | private fields read by `sigma`, `center`, `height` |
| `double region_rt_span_` (protected) | private field; `GaussTraceFitter::region_rt_span` (native getter) |
| `static const Size NUM_PARAMS_` (protected) | `GaussTraceFitter::NUM_PARAMS` |
| `void getOptimizedParameters_(const std::vector<double>&)` (protected) | `GaussTraceFitter::set_optimized_parameters([f64; 3])` (public) |
| `class GaussTraceFunctor` (protected): constructor, `operator()`, `df`, `m_data` | `GaussTraceFunctor<'a>`: `new(&MassTraces, weighted)`, `residuals`, `jacobian`, `inputs`, `values` (public) |
| `void setInitialParameters_(MassTraces&)` (protected) | `GaussTraceFitter::set_initial_parameters(&MassTraces) -> Result<()>` (public) |
| `void updateMembers_()` (protected) | `TraceFitter::set_parameters` |
| file-local `double pow2(double)` | private `pow2` |

The protected `setInitialParameters_`, `getOptimizedParameters_` and
`GaussTraceFunctor` are public here so that the start values, the stored
parameters and the residuals and Jacobians can be compared with the executed
C++ (the C2 oracle reaches them through a derived probe class).

## Preserved source conventions

`TraceFitter`:

- `optimize` refuses `values < x.len()` before evaluating anything, with the
  message "Skipping feature, we always expect N>=p".
- The solver is the Eigen `LevenbergMarquardt` with its defaults (`factor` 100,
  `ftol` = `xtol` = `sqrt(f64::EPSILON)`, `gtol` 0) and `maxfev` =
  `max_iteration`; the Jacobian is analytic and costs no function evaluation.
  `tests/trace_fitter.rs` checks the configuration against a direct call of
  `minimize` at every budget from 1 to 60.
- Every status after `ImproperInputParameters` is accepted, including
  `TooManyFunctionEvaluation`, and the parameters reached are copied back. A
  refused status gives "Could not fit the gaussian to the data: Error 0", with
  "gaussian" whichever model is fitted. `max_iteration <= 0` and an empty
  parameter vector reach that status without an evaluation, as Eigen's
  `minimizeInit` does. On every error the parameter vector is unchanged.
- `max_iteration` is signed, has no minimum and no maximum, and is read as an
  integer (through `from_param` only within `i32`; native difference 8);
  `weighted` is true exactly when it is the string `"true"`.
- `computeTheoretical` is `trace.theoretical_int * getValue(trace.peaks[k].rt)`.
- The span checks keep the subclasses' comparisons: `true` reports the
  violation, although the header documents the opposite.

`GaussTraceFitter`:

- The residual is `(baseline + theoretical_int * height * exp(c_fac * d^2) -
  intensity) * weight` with `c_fac = -0.5 / sigma^2` and `d = rt - x0`,
  evaluated left to right with the `f32` intensity promoted to `f64`; the
  Jacobian keeps `J0 = theo * e * w`, `J1 = theo * height * e * d * (1/sigma^2) * w`
  and the factor `0.125` in `J2 = 0.125 * theo * height * e * d^2 * (1/sigma^3) * w`,
  which is not the true derivative. Rows run trace by trace, then peak by peak.
- The weight is the theoretical intensity itself when `weighted`, otherwise 1.
- Start values: the intensity profile of `MassTraces::intensity_profile`; a
  centred five-point running mean over zero-padded totals, kept as a running
  sum in the source's order, skipped for `N <= 3`; the first strict maximum;
  `height = smoothed[max] - baseline`; `x0` at the maximum;
  `region_rt_span = last rt - first rt`; half-maximum walks with strict `>`
  that stop at the ends; `alpha = ((left + right) * 0.5) / height`;
  `sigma = 1.0` when `alpha >= 1`, otherwise `(delta_x * 0.5) / sqrt(-2 ln alpha)`,
  so a NaN `alpha` (for example an all-zero profile) gives a NaN `sigma`.
- After the fit `sigma` is stored as `|x[2]|`; height and centre unchanged.
- `getValue` computes `height * exp((-0.5 * d^2) / sigma^2)`, a different
  association from the residual's `c_fac * d^2`.
- The constants are the source's literals: FWHM `2.35482 * sigma`, area
  `2.506628 * height * sigma`, bounds `x0 -/+ 2.5 * sigma`,
  `checkMaximalRTSpan` `5.0 * sigma > max_rt_span * region_rt_span`,
  `checkMinimalRTSpan` `(upper - lower) < min_rt_span * 5.0 * sigma`.
- `getGnuplotFormula` writes `<name>(x)= <baseline> + <theo * height> *
  exp(-0.5*(x-<rt_shift + x0>)**2/(<sigma>)**2)` with the default stream
  precision 6.
- With one or two peaks the start values are set and then the fit is refused,
  so the fitter keeps the start values, as after the source's exception.

## Native differences

1. **`exp` and `log` come from the platform C library**, as in the source,
   through `f64::exp` and `f64::ln`, so fitted parameters are **not
   bit-identical across platforms**. glibc 2.39 (Linux x86-64) and Apple libm
   (macOS arm64) disagree in the last bit at 70 of the 37,080 distinct `exp`
   arguments the tests reach. The same commit therefore gives different fits on
   the two platforms the project tests on; for example `start.trailing_max`
   deviates from the oracle by 1.72e-3 on Linux and 1.49e-3 on macOS. Both meet
   every acceptance criterion. A fit is serial and repeats bit for bit on one
   platform.
   - The work package's notes proposed the `libm` crate for cross-machine
     reproducibility. It was measured and rejected: its `exp` differs from both
     platform libraries by one unit in the last place at 9% of the recorded
     residual points, which breaks the 1e-14 residual criterion and moves
     FeatureFinderCentroided_1 fits by up to 1.1e-9.
   - A correctly rounded `exp` and `log` meets every criterion in a
     table-lookup simulation and would be identical on every platform. It
     matches the oracle in fewer last bits.
   - The choice is the integrator's (see "`exp` and `log` across platforms").
     B5-EGH must follow it.
2. **Errors** are `Error::InvalidValue` whose message starts with
   `UnableToFit-FinalSet: `, the source exception's name; the part after it is
   the source's `what()`.
3. **Traces without any peak** are refused with the "N>=p" error before the
   fitter changes. The source reads the first element of the empty intensity
   profile first, which is undefined behaviour.
4. **`compute_theoretical`** returns `Error::InvalidValue` for an index past the
   peaks, where the source indexes without a check.
5. **Work ceilings.** `optimize` checks the solver's point and byte ceilings
   (`levenberg_marquardt::preflight_points`: 1,000,000 residuals, 64 MiB of
   Jacobian) before evaluating, and passes at most `MAX_RESIDUAL_WORK / values`
   (2^30 residual evaluations in total, at least one evaluation) as `maxfev`.
   When the ceiling is below the configured budget and the solver stops on
   it, the fit is refused with `Error::InvalidValue` and the parameters are
   left unchanged, where the source would continue or accept.
   - A fit that ends before the ceiling is exactly the uncapped fit.
   - So is a fit that ends at the ceiling's own evaluation with status 1, 2
     or 3, which the solver tests before the budget.
   - A fit whose uncapped run would end at exactly that evaluation with status
     4, 6, 7 or 8 is refused. Eigen and the transcription test the budget
     before `FtolTooSmall`, `XtolTooSmall` and `GtolTooSmall`, and before the
     next iteration's `CosinusTooSmall`. The unit tests
     `the_work_ceiling_refuses_only_fits_that_outlast_it` and
     `a_late_status_at_exactly_the_ceiling_is_refused` check both sides.

   Within the solver's point ceiling the work ceiling is at least 1,073
   evaluations. FeatureFinderCentroided_1's fits
   end after 30 to 131 evaluations over about 100 residuals, against a ceiling
   of about ten million.
6. **`Clone` copies `region_rt_span`.** The source copy constructor and
   assignment copy `height_`, `x0_` and `sigma_` but not `region_rt_span_`, so
   `checkMaximalRTSpan` on a copy reads an indeterminate value.
   `FeatureFinderAlgorithmPicked` never copies a fitter.
7. **Uninitialised members start at zero.** A new fitter's height, centre,
   sigma and region span are `0.0`; the source leaves them indeterminate until
   the first fit.
8. **`Param` mapping.** `from_param` fills defaults, validates types and valid
   strings as `Param::checkDefaults` does, and returns unknown keys as warnings;
   the typed record keeps only the two known values, where the source keeps the
   unknown key in `param_`. The wording of the validation errors is the
   `crate::param` module's, not the source's (for example "TraceFitter:
   Invalid string parameter value 'maybe' for parameter 'weighted' given!
   Valid values are: 'true,false'."); the executed cases in
   `tests/trace_fitter.rs` check the outcome, not the text.
   **`max_iteration` beyond `i32`:** `to_param` writes any `i64`, but
   `from_param` refuses a value outside `i32` ("parameter value cannot be
   converted to i32"), so the round trip fails there. The cause is
   `crate::param`'s restriction check (`src/param.rs`), which converts
   every integer entry to `i32`. The source's `ParamEntry::isValid`
   (`Param.cpp:122`) narrows to `int` without a check, and `updateMembers_`
   reads the full value. The executed SDK accepts and stores 2^31,
   3,000,000,000, -2^31 - 1 and `i64::MAX`
   (`../oracle/gauss-trace-fitter/param-range`). `FeatureFinderAlgorithmPicked`
   never reaches the range: the `UInt` it passes stays below 2^31, because the
   minimum of `fit:max_iterations` is checked on the `int`-narrowed value.
   `to_param_and_from_param_disagree_beyond_i32` records the current boundary.
9. **Borrowing.** `fit` borrows the traces immutably, and `area`, the span
   checks and `gnuplot_formula` take `&self`; the source declares them
   non-const, but no subclass writes.
10. **`gnuplot_formula`'s function name** is a Rust `char`; the source streams
    one C++ `char`. A character outside ASCII is written as its UTF-8 encoding.
    Numbers follow `format::file_info::text_format::ostream_g`, which rounds
    exact decimal ties as glibc does; Apple libc keeps trailing zeros on one
    class of ties (see that module).
11. **The status is available.** `optimize_with_status` returns the accepted
    `LmStatus`; the source discards it.

## Checked boundaries and evidence

| Evidence | Tier | Where |
| --- | --- | --- |
| The 16 `START_SECTION`s of `TraceFitter_test.cpp` | 3, and 4 for the compile-time form | `tests/trace_fitter.rs` (`section_*`); the `compile_fail` doctest on `TraceFitter` |
| The 17 `START_SECTION`s of `GaussTraceFitter_test.cpp`, literals and `TOLERANCE_RELATIVE(1.001)` | 3 | `tests/gauss_trace_fitter.rs` (`section_*`) |
| C2 class-test fits (theoretical intensities 0.8/0.2 and 0.4/0.6, weighted and unweighted): start values, residuals and Jacobians at five vectors, fit at 500, every query, budget sweep 1..500, evaluation path | 1 (library), adapted (replica status, `nfev`, `njev`, path) | `class_test_cases_replay_the_oracle`, `class_test_fits_match_the_c2_literals`, `class_test_budget_boundaries` |
| C2 degenerate inputs (`flat3`, `short3`, fewer peaks than parameters) | 1 / adapted | `degenerate_cases_replay_the_oracle` |
| The 25 Gaussian fits of FeatureFinderCentroided_1's seed loop, with the traces the loop handed the fitter | 1 / adapted | `feature_finder_centroided_1_seed_fits_replay_the_oracle` |
| `getDefaults`, eight `setParameters` cases, four fit failures and the state they leave, nine start-value boundaries, five gnuplot formatting cases | 1 | `tests/trace_fitter.rs`, `fit_failures_match_the_oracle`, `start_value_boundaries_match_the_oracle`, `gnuplot_number_formatting_matches_the_oracle` |
| `optimize` configuration, refusals, exhausted budget, huge budget, work ceiling, `initial_shape` loops, non-finite input | 4 | `tests/trace_fitter.rs`, unit tests in `trace_fitter.rs`, `non_finite_input_does_not_panic` |
| The review probe's 79 inputs: start vectors, region spans, residuals and Jacobians at the start vector | 1 (library), adapted for the replica's first evaluation point | `solver_gap_inputs_first_step_inputs_match_the_oracle` |
| The same inputs: fitted parameters, status, `nfev`, `njev`, residual path (**known gap**, see below) | 1 / adapted, recorded and reported, not asserted | `solver_gap_probe_reports_the_known_gap` (`#[ignore]`) |

**Oracle.** `../oracle/gauss-trace-fitter/run.sh` compiles `driver.cpp` against
the product SDK (Debug, core `4fdec46`; the traced sources are hash-identical to
`bc9cc12`) with `-ffp-contract=off`, runs it twice in the fixed oracle
environment and requires identical output, then extracts the Gaussian records of
the C2 class-level oracle (`../oracle/featurefinder-picked`, results `omp1` and
`omp4`, which must agree) with `extract_c2.py`. The four fixture files are its
output; their sha256 values are in the manifest. The fifth fixture,
`solver_gap.tsv`, comes from `../oracle/gauss-trace-fitter/solver-gap/run.sh`
(see "Known gap"). No C++ is in this repository.

**Comparison.** Start values, residuals, Jacobians and the queries (evaluated on
the oracle's fitted parameters) are compared within 1e-14 relative; fitted
parameters, sweep parameters and the evaluation path within 1e-9 relative;
statuses, `nfev`, `njev`, the evaluation order, booleans, strings and the budget
boundaries exactly. Measured against the macOS arm64 oracle on Linux x86-64
(IBMI dax, glibc 2.39) and natively on macOS arm64 with Apple libm. The two
libraries differ mainly in the bit-identical counts, through `exp` (native
difference 1).

Re-measured on 2026-09-15 for package B3b-LM-FIDELITY, which replaced the
solver's reduction arithmetic with Eigen's own kernels
(`docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1). Both columns are executed runs of
`cargo test --test gauss_trace_fitter -- --nocapture` on the merged tree: glibc
on dax, Apple libm on the development machine. Every deviation fell and every
bit-identical count rose:

| Replay | Values | Bit-identical, glibc | Bit-identical, Apple libm | Largest relative deviation (both unless noted) |
| --- | ---: | ---: | ---: | --- |
| class-test cases | 10,461 | 8,337 | 8,338 | start values, residuals, Jacobians, queries 0; fit 1.44e-12; sweep 1.81e-12; path 9.08e-12 |
| FeatureFinderCentroided_1 seeds | 19,572 | 17,679 | 17,930 | start values, residuals, Jacobians, queries 0; fit and sweep 1.83e-10 glibc, 1.84e-10 Apple libm |
| degenerate inputs | 24 | 24 | 24 | 0 |
| fit failures | 32 | 32 | 32 | 0 |
| gnuplot case | 3 | 3 | 3 | 0 |

Before B3b the same replays gave: class tests 8,255 / 8,255 bit-identical with
fit 2.5e-12, sweep 3.7e-12 and path 1.8e-11; seeds 16,699 / 16,696 with fit and
sweep 6.4e-10; degenerate inputs 22 / 22 with fit 3.0e-15. The start-value
boundary cases moved in both directions and are treated separately below.

On these fixtures, every status agrees at every budget checked: 1 to 500 for
the class-test cases, and for the seeds every budget up to two past natural
termination plus 100, 131, 250, 499 and 500. `nfev` and `njev` agree at every
class-test budget and, for the seeds, at 500; the sequence of residual and
Jacobian evaluations of the four class-test fits at 500 is the same, with
parameters within 1e-9. **This agreement is fixture-specific.** It is not a
general property of the port: on other inputs the fits, statuses and `nfev`
differ from the executed solver (see "Known gap: solver fidelity beyond the
fixtures"). The
budget boundaries are the
oracle's: the class-test results equal the result at 500 from `max_iteration`
47 (0.8/0.2, unweighted and weighted) and 31 (weighted 0.4/0.6) and differ at
every smaller budget; each seed's boundary is met exactly, 29 to 46 for the 24
seeds that become features and 131 for seed 24, whose fit the algorithm then
rejects ("Fitted model is bigger than 'max_rt_span'"). The work package's literals hold within 1e-9: height
9.9995822835781674 and sigma 1.5000586589470672 unweighted, height
6.0828584222534152 weighted 0.4/0.6.

**Ill-conditioned boundary cases.** Four of this package's own start-value
cases, not work-package targets, carry a measured bound instead of the 1e-9 fit
tolerance. Two of them still need it after B3b-LM-FIDELITY.

What was executed:

- Their start values and intensity profiles are bit-identical to the
  oracle's.
- The driver did not record the C++ functor for these cases, so their
  residuals and Jacobians were not compared.
- Before B3b the gap was the same with glibc's `exp` and `log`, with Apple
  libm's (the oracle platform's library: the macOS review run, and the
  substitution run on Linux) and with correctly rounded ones, except for
  `trailing_max`. After B3b the two remaining cases agree to twelve digits
  between glibc and Apple libm and `trailing_max` agrees exactly; only
  `merged_profile` splits, and only by being bit-identical on macOS where glibc
  leaves 2.69e-11.

So the gap does not come from `exp` or `log`. These four cases are not
special: they are instances of the general solver gap described in the next
section, which the review's 79 further inputs exposed.

The bounds were **re-measured on 2026-09-15** for package B3b-LM-FIDELITY,
which replaced the solver's reduction arithmetic with Eigen's own kernels
(`docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1). Both columns below are that
package's own executed runs of
`start_value_boundaries_match_the_oracle -- --nocapture`: Linux x86-64 with
glibc 2.39 on dax, and macOS arm64 with Apple libm on the development machine.

| Case | Deviation before B3b: glibc / Apple libm / correctly rounded | Deviation after B3b: glibc / Apple libm | Asserted bound |
| --- | --- | --- | --- |
| `start.n4_boundary` (four points, maximum second) | 1.40e-9 / 1.40e-9 / 1.40e-9 | bit-identical / bit-identical | 1e-8 |
| `start.merged_profile` (two traces with interleaved retention times) | 4.41e-9 / 4.41e-9 / 4.41e-9 | 2.69e-11 / bit-identical | 1e-8 |
| `start.leading_max` (maximum at the first retention time, sigma 0.25) | 2.93e-4 / 2.93e-4 / 2.93e-4 | **2.06e-3 / 2.06e-3** | **1e-2** (was 1e-3) |
| `start.trailing_max` (maximum at the last retention time) | 1.72e-3 / 1.49e-3 / 1.49e-3 | 1.72e-3 / 1.72e-3 | 1e-2 |

The expected values are the oracle's. The bounds are not: each is the next
power of ten above the port's own largest measured deviation. The first two
now sit under the package's own 1e-9 fit tolerance on both measured platforms;
their 1e-8 bounds are kept as headroom for the targets neither package could
measure (the `cross-platform` CI job also builds Windows), not because the
deviation needs them.

**Why `leading_max` rose, and why the bound was raised.** This is the one
number in this package that the new solver moves away from its oracle, and it
is a consequence of a decision taken by the user on 2026-09-15 and recorded in
`docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1: the port matches Eigen as the
**Linux x86_64 Release** build compiles it, on every target. Eigen's arm64
kernels fuse their `pmadd` lanes (`vfmaq_f64` under `EIGEN_VECTORIZE_FMA`,
which `__ARM_FEATURE_FMA` always sets) and the x86_64 `-mssse3` kernels do not,
so one unfused Rust path cannot reproduce both C++ builds. This package's
fixtures were generated by the **macOS arm64** product SDK, which is the fused
build. On the traced-path oracle the committed solver keeps 141 of 141
evaluation paths and final parameters against Linux x86_64 Release and 21 of
141 against the macOS SDK.

`start.leading_max`'s fit is ill-conditioned enough to make that visible:
2.0614474501769614e-3 relative on Linux and 2.0614474501789866e-3 on macOS,
against 2.93e-4 before, so the previous 1e-3 bound would fail. The bound moves
to the next power of ten, 1e-2, and is the same bound `trailing_max` already
carried. No oracle value and no other tolerance in this package changed.

**Everything else in this package got closer or stayed inside its bound**,
measured in the same runs. `start.trailing_max` is the only other movement: on
macOS arm64 it went from 1.49e-3 to the 1.72e-3 Linux already measured, well
inside its unchanged 1e-2 bound.

| Quantity | Before B3b | After B3b (dax) |
| --- | --- | --- |
| Class-test fitted parameters | 2.51e-12 | 1.44e-12 |
| Class-test evaluation path (`path.x`) | 1.79e-11 | 9.08e-12 |
| FeatureFinderCentroided_1 seed fits | 6.36e-10 | 1.83e-10 |
| FeatureFinderCentroided_1 seed values bit-identical | - | 17,679 of 19,572 |
| `start.n4_boundary`, `start.merged_profile` | 1.40e-9, 4.41e-9 | bit-identical, 2.69e-11 |

**Known gap: solver fidelity beyond the fixtures.** `GaussTraceFitter::fit`
is not bit-faithful to the **macOS arm64** library in general. The class-test
and FeatureFinderCentroided_1 results above (fits within 1.44e-12 and 1.83e-10,
equal statuses, `nfev` and `njev`) are fixture-specific.

The root cause this section left open was found and fixed in package
B3b-LM-FIDELITY: the transcription summed left to right where Eigen accumulates
in SIMD lanes. `docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1 has the traced
first-divergence evidence and `tests/lm_eigen_path_differential.rs` the tier-1
assertion. What remains is not a defect but a platform split that the C++ has
too: Eigen's arm64 kernels fuse their `pmadd` lanes and the x86_64 `-mssse3`
kernels do not, and on 2026-09-15 the user decided the port matches the **Linux
x86_64 Release** build on every target. This section's oracle was generated by
the macOS arm64 SDK, the fused build, so a residual gap against it is expected
and is the price of that decision. Against Linux x86_64 Release the solver is
bit-exact on all 141 traced fits.

What was executed:

- The round-2 review of this package (commit `e6a5fec`) generated 79 inputs:
  12 edge cases, 7 budget cases (`max_iteration` 1 to 13) and 60 seeded random
  Gaussian sets (1 to 4 traces, 3 to 40 peaks, 0 to 30% noise, weighted or
  unweighted, budget 1 to 60 or 500).
- On macOS arm64, the oracle platform with Apple libm, it ran them through
  the product-SDK `GaussTraceFitter::fit` (tier 1; two runs, byte-identical
  output). It also ran them through its Eigen 5.0.1 replica of `optimize_` on
  the library's own `GaussTraceFunctor` and start values (adapted: status,
  `nfev`, `njev`, the raw final vector and every residual-evaluation point).
  The replica reproduces the library's final parameters bit for bit in 79 of
  79 cases, so its path is the library's.
- The review's drivers are promoted unchanged to
  `../oracle/gauss-trace-fitter/solver-gap` (`gen_cases.py`, `replica.cpp`,
  and `review_probe.cpp` as `library_probe.cpp`). This round added
  `functor_probe.cpp`: the library's `GaussTraceFunctor` residuals and
  Jacobian at the library's start vector, the inputs of Eigen's first step
  (tier 1). `run.sh` builds all three like the other drivers
  (`-O0 -ffp-contract=off`, product SDK, Eigen 5.0.1), runs each twice with
  byte-identical output, and `to_fixture.py` joins the results. Rerun on
  2026-09-14, it reproduced the review's library and replica outputs byte for
  byte. The fixture `tests/data/gauss_trace_fitter/solver_gap.tsv` holds the
  inputs and the executed results, and nothing computed by Rust.
- `solver_gap_inputs_first_step_inputs_match_the_oracle` asserts the start
  vectors, the region spans, and the residuals and Jacobians at the start
  vector within 1e-14. All 14,065 values are bit-identical on both platforms.
- `solver_gap_probe_reports_the_known_gap` is `#[ignore]`d. Run with
  `--ignored --nocapture`, it prints the port's deviations case by case and
  asserts nothing about them. The printed reports of the two runs summarised
  below are kept in `../oracle/gauss-trace-fitter/solver-gap/logs`.

Re-measured on 2026-09-15 with the B3b-LM-FIDELITY solver: the macOS arm64
column is a native run on the development machine, the Linux column a run on
dax (glibc 2.39), both of
`cargo test --test gauss_trace_fitter solver_gap_probe_reports_the_known_gap
-- --ignored --nocapture`. The "before" column is the pre-B3b measurement this
section carried, which was the same on both platforms except where noted. These
are measurements of the current port against a macOS-generated oracle, not
expectations:

| Measure | Before B3b | macOS arm64 | Linux x86-64 |
| --- | --- | --- | --- |
| Start vector (evaluation 0) bit-identical | 79 of 79 | 79 of 79 | 79 of 79 |
| Residual path departs from Eigen's | 72 of 79: 58 at evaluation 1 | 64 of 79: 30 at evaluation 1, the rest by evaluation 43 | 65 of 79: 30 at evaluation 1, the rest by evaluation 73 |
| Largest coordinate's relative deviation at the first differing evaluation | 1.1e-16 to 8.3e-14 | 1.1e-16 to 2.6e-14 | the same |
| Final height, centre, sigma and span bit-identical | 124 / 123 of 316 | 151 of 316 | 145 of 316 |
| Fits beyond 1e-9, natural termination | 12 | 12, by 1.04e-9 to 1.16e-4 | 12, by 1.06e-9 to 9.52e-5 |
| Fits beyond 1e-9, budget of 500 exhausted | 9 | 8, by 2.34e-7 to 2.04e-2 | 9, by 1.64e-7 to 3.72e-2 |
| Fits stopped by a budget below 500 (20) | within 4.9e-12 | within 8.48e-13 | within 8.48e-13 |
| Status differs | 3 / 2 | 1 (`random_52`) | 1 (`random_52`) |
| `nfev` differs at natural termination | 9 | 7 | 6 |
| `njev` differs | 17 | 14 | 13 |

The fits beyond 1e-9 (relative deviation of height, centre or sigma from the
library). `large_intensity`, which was the worst natural termination before at
6.6e-5 with `nfev` 188 against Eigen's 249, is now within 1e-9 and gone from
this table; `random_34` entered it at 1.14e-9:

| Case | Eigen status/`nfev` | macOS arm64 | Linux x86-64 |
| --- | --- | --- | --- |
| `negative_intensity` | 1/132 | 3.62e-7, 1/**131** | 9.45e-8, 1/132 |
| `unsorted_two_traces` | 1/118 | 8.74e-7, 1/118 | 8.74e-7, 1/118 |
| `random_7` | 1/86 | 1.04e-9, 1/86 | 1.06e-9, 1/86 |
| `random_17` | 2/378 | 1.16e-4, 2/**388** | 3.73e-5, 2/**389** |
| `random_21` | 1/112 | 8.24e-5, 1/**114** | 8.24e-5, 1/**114** |
| `random_28` | 2/126 | 8.56e-5, 2/**129** | 6.72e-5, 2/**129** |
| `random_33` | 1/97 | 6.41e-9, 1/97 | 6.41e-9, 1/97 |
| `random_34` | 1/75 | 1.14e-9, 1/75 | 1.14e-9, 1/75 |
| `random_39` | 1/314 | 7.70e-6, 1/**318** | 7.70e-6, 1/**318** |
| `random_48` | 3/52 | 5.29e-9, 3/52 | 5.29e-9, 3/52 |
| `random_52` | 2/157 | 9.52e-5, **1**/**146** | 9.52e-5, **1**/**146** |
| `random_55` | 1/111 | 7.19e-5, 1/**131** | 7.19e-5, 1/**131** |
| `weighted_theo_negative` | 5/500 | within 1e-9 | 1.01e-2 |
| `random_4` | 5/500 | 4.23e-4 | 6.30e-5 |
| `random_8` | 5/500 | 1.29e-3 | 3.62e-4 |
| `random_10` | 5/500 | 2.26e-4 | 8.55e-4 |
| `random_12` | 5/500 | 3.08e-3 | 3.33e-3 |
| `random_23` | 5/500 | 3.98e-3 | 2.43e-3 |
| `random_24` | 5/500 | 2.34e-7 | 1.64e-7 |
| `random_32` | 5/500 | 2.20e-3 | 1.65e-3 |
| `random_58` | 5/500 | 2.04e-2 | 3.72e-2 |

What this establishes:

- On all 79 inputs the start vector and the residuals and Jacobian there are
  bit-identical to the library's, on both platforms, before and after B3b.
- **The root cause was in `src/math/fitters/levenberg_marquardt.rs` and is
  fixed:** the transcription summed left to right where Eigen accumulates in
  SIMD lanes. `docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1 records the traced
  first divergence of 141 fits, the kernels that now reproduce it, and the
  tier-1 assertion in `tests/lm_eigen_path_differential.rs`.
- **What remains against this oracle is the platform split, not a defect.**
  Eigen's arm64 kernels fuse their `pmadd` lanes; the x86_64 `-mssse3` kernels
  do not; this fixture was generated by the fused macOS arm64 build, and the
  user decided on 2026-09-15 that the port matches the unfused Linux x86_64
  Release build everywhere. Against a Linux x86_64 Release oracle the solver
  reproduces all 141 traced fits bit for bit, evaluation path included.
- The two ill-conditioned driver cases above are instances of this split, not
  its extent. The 21 fits in the table are the cases where an ill-conditioned
  or non-converged problem amplifies a last-bit difference; the two largest
  (`random_58`, `weighted_theo_negative`) are both status-5 fits that
  exhausted their 500-evaluation budget and never converged.
- `exp` and `log` are second-order here, as before.
- Asserting this report is not the right follow-up while the fixture is
  macOS-generated. The asserted replay lives in
  `tests/lm_eigen_path_differential.rs`, where the same 79 inputs are compared
  against a **Linux x86_64 Release** trace of Eigen and are bit-identical.

**`exp` and `log` across platforms.** The evidence for native difference 1:

- **Recorded points.** At the 840 class-test and 3,132 seed residual points
  the oracle recorded, glibc's `exp` (dax) and Apple libm's `exp` (macOS arm64,
  review run) reproduce every C++ residual and Jacobian entry bit for bit.
- **The `libm` crate.** Its `exp` differed from the platform libraries in 22
  and 340 of those evaluations, by one unit in the last place. With it, 21 and
  204 residuals were not bit-identical, and the largest residual deviation was
  9.1e-12 relative (a residual near convergence). The
  FeatureFinderCentroided_1 fits reached 1.13e-9 relative in sigma for 3 of 25
  seeds, and the class-test `checkMaximalRTSpan` knife edge
  (`5 * sigma / span + 1e-14`) flipped.
- **Every argument the tests reach**
  (`../oracle/gauss-trace-fitter/exp-simulation`, `manifest.json` sha256
  `09fb44d732ca9f318d0df8c25f1be5f1e50f85d13429c81fe9646699259a9da7`). A
  scratch checkout of `46be7d2` recorded every distinct `exp` and `log`
  argument of `tests/gauss_trace_fitter.rs` on dax. It then
  replaced both functions with a lookup table and repeated the tests until no
  argument was missing, once with Apple libm's values (CPython `math` on macOS
  arm64) and once with correctly rounded values (80-digit `decimal`). Over the
  37,080 distinct `exp` and 25 `log` arguments of the three runs:

  | `exp` values | Disagreements |
  | --- | ---: |
  | glibc 2.39 vs Apple libm | 70 |
  | glibc vs correctly rounded | 31 |
  | Apple libm vs correctly rounded | 69 |
  | misrounded identically by both libraries (error 0.5000-0.502 ulp) | 15 |
  | `log`, any pair | 0 |

  The three substitutions give these replays, all 32 tests passing:

  | `exp` and `log` | Class tests bit-identical | Seeds bit-identical | `tie4_smoothed` | `flat_alpha_ge_1` | `trailing_max` |
  | --- | ---: | ---: | --- | --- | --- |
  | glibc (dax) | 8,255 | 16,699 | 3.74e-13 | 2.95e-11 | 1.72e-3 |
  | Apple libm substituted on dax | 8,255 | 16,696 | 0 | 5.95e-15 | 1.49e-3 |
  | correctly rounded, substituted on dax | 6,355 | 16,658 | 3.74e-13 | 2.95e-11 | 1.49e-3 |

  With Apple libm substituted, the Linux run reproduces the macOS review run
  (16,696 seed values; `tie4_smoothed` 0, `flat_alpha_ge_1` 5.95e-15,
  `trailing_max` 1.49e-3), so the C library is the whole cross-platform
  difference.
- **Correctly rounded `exp` and `log`.** Residuals and Jacobians still deviate
  0 at every recorded point, and every other acceptance criterion holds. The
  bit-identical counts fall because neither platform library is correctly
  rounded.
- **Decision.** Keeping the platform library matches the oracle platform's
  last bits; a correctly rounded implementation gives the same bits on every
  platform. Choosing is the integrator's decision (a crate such as a
  CORE-MATH port would need a `Cargo.toml` change). B5-EGH must make the same
  choice.

**Not covered.** A negative `max_iteration` through the tool path cannot occur
(`fit:max_iterations` has minimum 1). The C++ copy constructor's uninitialised
`region_rt_span_` is undefined behaviour and was not executed.

## C++ issue candidates

1. **`GaussTraceFitter` copy leaves `region_rt_span_` uninitialised**
   (`GaussTraceFitter.cpp:25-46`): the copy constructor and assignment copy
   three of the four model members, so `checkMaximalRTSpan` on a copy reads an
   indeterminate value. Source review; not executed.
2. **The sigma column of the Jacobian carries `0.125`**
   (`GaussTraceFitter.cpp:199`), not the derivative's `1`. The fit still
   converges on the tested data, but the Levenberg-Marquardt path, and with it
   the budget boundaries, depend on the wrong derivative. Source review; the
   executed Jacobian confirms the factor.
3. **An all-zero intensity profile "fits" with a NaN sigma**
   (`GaussTraceFitter.cpp:283-291` and `TraceFitter.cpp:127-130`): height 0
   gives `alpha = 0/0`, `sigma = NaN`, NaN residuals, and Eigen stops with
   `CosinusTooSmall`, which `optimize_` accepts. Executed (`start.all_zero`).
4. **Empty traces are undefined behaviour** in `setInitialParameters_`
   (`smoothed[max_index]` on an empty vector), and `computeTheoretical` indexes
   `trace.peaks[k]` unchecked. Source review.
5. **Documentation.** `TraceFitter.h` (lines 170-191) documents
   `checkMinimalRTSpan` and `checkMaximalRTSpan` with the opposite meaning of
   both implementations and of
   `FeatureFinderAlgorithmPicked`'s use, and calls `max_iteration` a number of
   iterations although it is passed as `maxfev`. `optimize_` reports "Could not
   fit the gaussian" for the EGH model too.

## Ledger notes

- `TraceFitter.h`: every public and protected member is mapped above; proposed
  `complete` with `tests/trace_fitter.rs` and this document.
- `GaussTraceFitter.h`: every member is mapped above; proposed `complete` with
  `tests/gauss_trace_fitter.rs` and this document.
- Known gap for the ledger: `GaussTraceFitter::fit` is bit-faithful to Eigen
  as the Linux x86_64 Release build compiles it, and not to the macOS arm64
  product SDK that generated these fixtures (see "Known gap: solver fidelity
  beyond the fixtures"). The transcription defect that lane B3b owned in
  `src/math/fitters/levenberg_marquardt.rs` is fixed; what remains is the
  platform split Eigen has, resolved by the user's decision of 2026-09-15 in
  favour of the Linux build (`docs/DISTRIBUTION_FITTERS_SUPPORT.md` §1). It
  applies to every consumer of that solver, including B5-EGH and the feature
  finder's fits, and it is not limited to the ill-conditioned cases.
- New module edges `analysis -> math` and `analysis -> param`, both acyclic;
  `tools/check_module_cycles.py` reports them as not yet recorded.
- CI: `cargo test --locked --no-default-features --test trace_fitter --test gauss_trace_fitter`
  in the `minimum-rust` job.
