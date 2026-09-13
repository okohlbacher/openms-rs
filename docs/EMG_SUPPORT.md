# Exponentially modified Gaussian fitting

`analysis::emg::EmgGradientDescent` ports the scientific operations of OpenMS
`MATH/MISC/EmgGradientDescent` from revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It fits native spectra and
chromatograms, estimates parameters from `f64` slices, and evaluates supplied
parameters with optional extrapolation. It does not use C++ or an external
optimizer. `libm` 0.2.16 supplies only the complementary error function; the
remaining mathematical functions use Rust's standard library.

**Target revision.** The header, the implementation and the class test are
byte-identical at the current target `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
and at the original pin: `EmgGradientDescent.h` is
`91586579fd7e311fe56479a1397b1870ffd0c84aa634b135ac4fd27129dc9303`,
`EmgGradientDescent.cpp` is
`584a23ffa22efdba082060f6160c6a4b2c85a697c7510e32fc488db182bc613e` and
`EmgGradientDescent_test.cpp` is
`62832ae3a40e258f4cb2480f6b2c5899310bfcf6c51ed11654746d09a61b8019` in both
checkouts, so the port carries forward unchanged. `tests/data/emg_provenance.json`
records that as an empty `target_verification.source_changes`.

Rust files covering the header: `src/analysis/emg.rs` (one module). Tests:
`tests/emg.rs`, `tests/emg_reference.rs`, `tests/emg_integration.rs`,
`tests/emg_sections.rs` and five private unit tests in the module.

## API mapping

Every member of the header appears below, public, protected and private. The
C++ private and protected members are reachable only through the header's own
`EmgGradientDescent_friend` shim, which exists purely so the class test can
call them; the Rust counterparts are private functions covered by unit tests in
the same file, which is the same arrangement without the shim.

| C++ member | Visibility | Rust counterpart | Notes |
|---|---|---|---|
| `EmgGradientDescent()` | public | `EmgGradientDescent::default()` | Constructs with the source's own defaults instead of routing them through `Param`. |
| `~EmgGradientDescent() override = default` | public | derived drop | The value owns nothing needing a destructor. |
| `void getDefaultParameters(Param& params)` | public | `impl Default for EmgGradientDescent` | The three `Param` entries become two typed fields plus a dropped one; see "native differences". |
| `friend class EmgGradientDescent_friend` | public | not ported: a test-only access shim. Rust unit tests live in the module and need no friendship. |
| `template<PeakContainerT> void fitEMGPeakModel(const PeakContainerT&, PeakContainerT&, double left_pos = 0.0, double right_pos = 0.0) const` | public | `fit_spectrum(&MSSpectrum, Option<f64>, Option<f64>) -> Result<EmgSpectrumFit>` and `fit_chromatogram(&MSChromatogram, …) -> Result<EmgChromatogramFit>` | The two explicit template instantiations become two named methods. The out-parameter becomes the return; `0.0` as "no bound" becomes `None`, so a literal zero bound is expressible. |
| `UInt estimateEmgParameters(const vector<double>& xs, const vector<double>& ys, double& best_h, double& best_mu, double& best_sigma, double& best_tau) const` | public | `estimate_parameters(&[f64], &[f64]) -> Result<EmgEstimate>` | Four out-parameters and the `UInt` iteration count all return in `EmgEstimate`. |
| `void applyEstimatedParameters(const vector<double>& xs, double h, double mu, double sigma, double tau, vector<double>& out_xs, vector<double>& out_ys) const` | public | `apply_parameters(&[f64], EmgParameters) -> Result<EmgCurve>` | Four scalars become `EmgParameters`; two out-parameters become `EmgCurve`. |
| `void updateMembers_() override` | protected | not ported: there is no `Param` indirection to react to. Fields are read directly. |
| `void extractTrainingSet(const vector<double>& xs, const vector<double>& ys, vector<double>& TrX, vector<double>& TrY) const` | protected | private `training_set(&[f64], &[f64]) -> Result<(Vec<f64>, Vec<f64>)>` | Collection order preserved exactly; `Exception::SizeUnderflow` becomes `Error::InvalidValue`. |
| `double computeMuMaxDistance(const vector<double>& xs) const` | protected | private `mu_max_distance(&[f64]) -> f64` | `minmax_element` reproduced, including the `0.0` an empty container yields. Called with the training positions, as in the source. |
| `double computeInitialMean(const vector<double>& xs, const vector<double>& ys) const` | protected | private `initial_mean(&[f64], &[f64]) -> Result<f64>` | The six percentage levels and the running left/right pointers are transcribed; `Exception::SizeUnderflow` becomes an error. |
| `void iRpropPlus(double prev_diff_E_param, double& diff_E_param, double& param_lr, double& param_update, double& param, double current_E, double previous_E) const` | private | private `irprop_plus(f64, &mut f64, &mut f64, &mut f64, &mut f64, f64, f64)` | All four in/out references kept, because the caller owns four independent parameter states. |
| `double Loss_function(const vector<double>&, const vector<double>&, double h, double mu, double sigma, double tau) const` | private | private `loss(&[f64], &[f64], EmgParameters) -> Result<f64>` | Per-point division by `xs.size()` before summation, as in the source. |
| `double E_wrt_h(…)` | private | private `gradient_h(&[f64], &[f64], EmgParameters) -> Result<f64>` | Three branches transcribed literally. |
| `double E_wrt_mu(…)` | private | private `gradient_mu(…)` | As above. |
| `double E_wrt_sigma(…)` | private | private `gradient_sigma(…)` | As above. |
| `double E_wrt_tau(…)` | private | private `gradient_tau(…)` | As above. Note the source redeclares a local `PI` in this one function only; it holds the same value. |
| `double compute_z(double x, double mu, double sigma, double tau) const` | private | private `z(f64, EmgParameters) -> f64` | `(1/sqrt(2)) * (sigma/tau - (x-mu)/sigma)`, unchanged. |
| `double emg_point(double x, double h, double mu, double sigma, double tau) const` | private | private `model(f64, EmgParameters) -> Result<f64>` | Three branches transcribed; a non-finite result becomes an error instead of propagating. |
| `const double PI = OpenMS::Constants::PI` | private | `std::f64::consts::PI` | Same value; `Constants::PI` is itself the standard constant. |
| `UInt print_debug_` | private | not ported: no terminal output. The information it prints (iteration counts, training size, best loss and iterate) is returned in `EmgEstimate`. |
| `UInt max_gd_iter_` | private | `EmgGradientDescent::max_iterations` | |
| `bool compute_additional_points_` | private | `EmgGradientDescent::compute_additional_points` | |
| `class EmgGradientDescent_friend` and its nine forwarding methods | public (test shim) | not ported: private Rust functions are directly callable from the module's own `#[cfg(test)]` block. |

Native additions with no source counterpart: `EmgParameters`, `EmgEstimate`,
`EmgCurve`, `EmgSpectrumFit`, `EmgChromatogramFit`, `EmgParameters::validate`,
`EmgGradientDescent::validate`, and the `max_points` / `max_evaluations`
ceilings.

## Native API

```rust
use openms::analysis::emg::{EmgGradientDescent, EmgParameters};

let fitter = EmgGradientDescent {
    compute_additional_points: false,
    ..Default::default()
};
let parameters = EmgParameters { h: 100.0, mu: 10.0, sigma: 1.0, tau: 1.0 };
let curve = fitter.apply_parameters(&[9.0, 10.0, 11.0], parameters)?;
assert_eq!(curve.positions.len(), 3);
# Ok::<(), openms::Error>(())
```

`estimate_parameters(xs, ys)` returns `EmgEstimate`, containing the best
`EmgParameters { h, mu, sigma, tau }`, training mean squared error, training
point count, evaluated iteration count, one-based best iteration, scalar
evaluation count and a convergence flag. Exhausting `max_iterations` returns
the best finite result with `converged == false`; convergence means only the
source loss-history criterion was met, not evidence of a global optimum.

`apply_parameters(xs, parameters)` returns `EmgCurve` with `f64` positions and
intensities. `fit_spectrum(input, left, right)` and `fit_chromatogram(...)`
return a fitted native container, its estimate, and names of omitted data
arrays. Optional bounds are inclusive; `None` is unbounded and `Some(0.0)`
means exactly zero. Record metadata and container representation settings are
preserved. Fitted container intensities are checked before rounding to `f32`.

The C++ method creates a four-value `emg_parameters` FloatDataArray regardless
of peak count. Native parameters remain in the typed result because that array
would violate the kernel's alignment requirement. Input float, integer and
string data arrays are omitted and reported: no source rule defines their
values on changed or extrapolated samples. Inputs are never modified, including
when a later calculation fails.

## Preserved source algorithm

The defaults are 100,000 iterations and additional points enabled. The source
debug-print option is represented by returned diagnostics rather than logging
controls. Positions retain their input units; the fitter does not convert
minutes to seconds, center coordinates, or normalize intensities.

The initial amplitude is the maximum intensity. The initial mean averages six
sampled midpoint estimates at 60%, 65%, 70%, 75%, 80% and 85% of that maximum.
Sigma starts at one percent of this **absolute mean**, and tau at twice sigma.
The training set combines outer samples below 80% with selected inner slopes.
Its collection order is retained: left samples, right samples in reverse, then
selected inner left and right samples. It is not sorted before accumulation.

The source's three EMG expressions and all four analytical loss gradients are
transcribed with their original grouping, floating-point powers and per-point
normalization. iRprop+ starts each rate at 0.0125, uses factors 1.2 and 0.5,
caps rates at 2000, and undoes a previous update on a sign change when loss
increased. An exactly zero gradient still decreases the parameter by its
learning rate, as explicitly tested in OpenMS.

After simultaneous updates, amplitude cannot fall below the original maximum,
mean remains within 35% of the training span around its initial value, sigma is
clamped to `[0.0001, 20]`, and tau to `[sigma, 15*sigma]`. Every 50 iterations,
the current loss enters a ten-element ring initially filled with zeros. The
population standard deviation below 1 stops fitting. Parameters from the first
strictly best loss are retained. Native diagnostics count evaluated iterations;
they omit the C++ extra increment after exhausting the iteration limit.

Extrapolation uses the arithmetic mean of consecutive input spacings. It
extends only the side with the higher modeled endpoint, until intensity reaches
the other endpoint's intensity or 0.001, or before its distance from the sampled
apex exceeds three times the opposite side's distance. A generated point which
crosses an intensity threshold remains included. Equal modeled endpoints do
not cause extrapolation. Left extension preserves the same coordinates and
order without repeated insertion at the front of a vector.

## Checked differences and limits

- Every coordinate and intensity must be finite; coordinates must be strictly
  increasing. Estimation requires at least two points and equal slice lengths.
  Signed intensities are retained. The initial mean must be positive so the
  source initialization gives positive sigma and tau. Supplied sigma and tau
  must be positive, while amplitude and mean may be any finite values.
- Applying supplied parameters without extrapolation accepts empty and
  single-point input. Extrapolation requires at least two points.
- Zero iterations and zero resource limits are errors. Nonfinite model values,
  gradients, loss, updates, convergence statistics, generated coordinates, and
  unrepresentable `f32` intensities return errors. This includes failures after
  an earlier finite best: C++ can instead stop and return the previous best.
  No partially fitted result or `DBL_MAX` sentinel is returned.
- The source's middle model branch directly evaluates `exp(z*z) * erfc(z)` up
  to `z == 6.71e7`. It can overflow at much smaller positive values near 27.
  This port detects failure; it does not silently replace the model or its
  gradients with a scaled complementary error function or a new approximation.
- `max_points` defaults to 1,000,000 and bounds the whole input before container
  validation/allocation, as well as the generated output. `max_evaluations`
  defaults to 100,000,000 scalar model/gradient evaluations per call. Each
  optimizer iteration charges five evaluations per training point; application
  charges every original and generated point. A container fit shares this
  budget across estimation and application. Estimate diagnostics report only
  the estimation portion. Limits are checked before the corresponding work.

For `n` input points, `t <= n` training points, `i` evaluated iterations and
`k` generated points, numerical work is `O(n + i*t + k)` and working storage
is `O(n + k)`, plus copied record metadata and input annotations. Validation
and training extraction are linear and separately bounded by `max_points`.
Parameter estimation is not invariant to a change of coordinate units; the
source's minute and second examples deliberately use different initial values.

## Integration and evidence

`analysis::peak_integrator::PeakIntegrator` accepts optional EMG fitting.
Source preprocessing selects the requested interval, fits and extrapolates,
then replaces the interval with the fitted span for integration, background
estimation and shape metrics. Caller-supplied apex and height remain unchanged.
The ordinary default continues to integrate observed samples.

### Class-test section coverage

`EmgGradientDescent_test.cpp` has **13** `START_SECTION` blocks and all 13 are
mapped. Each row cites the Rust test and one concrete value it reproduces.

| Section | Rust test | Cited value |
|---|---|---|
| `EmgGradientDescent()` | `emg_sections::default_construction_reproduces_the_source_parameter_defaults` | `max_iterations == 100_000` |
| `~EmgGradientDescent()` | same test | the cloned value still reads `100_000` after the original is dropped |
| `getParameters()` | same test | `compute_additional_points == true` |
| `fitEMGPeakModel(MSChromatogram)` | `emg_reference::seven_literal_fits_match_source_parameters_counts_and_container_overloads` | `cutoff_min` fits `h == 3791.07` and produces 28 output points |
| `fitEMGPeakModel(MSSpectrum)` | same test, which runs both container overloads per trace | `saturated_min` produces 87 points on the spectrum overload too |
| `Loss_function(...)` | `emg_reference::seven_literal_fits_match_source_parameters_counts_and_container_overloads`, whose last assertion per trace is the source's own loss convention | `cutoff_min` full-trace loss `651.824632922326` |
| `extractTrainingSet(...)` | `emg_reference::source_initial_mean_and_training_count_are_observable_before_updates` and `emg_reference::ordered_training_selection_matches_independent_source_translation` | `saturated_min` selects 77 of its 83 points |
| `computeMuMaxDistance(...)` | private `emg::tests::source_mu_max_distance_uses_the_unsorted_span_and_tolerates_an_empty_set`, plus `emg_sections::the_fitted_mean_stays_within_thirty_five_percent_of_the_training_span` | `mu_max_distance([3,2,4,2,4,5,7,9,3]) == 2.45` and `mu_max_distance([]) == 0.0` |
| `computeInitialMean(...)` | `emg_reference::source_initial_mean_and_training_count_are_observable_before_updates` | `glutamate` initial mean `2.69743333333333` |
| `iRpropPlus(...)` | private `emg::tests::source_irprop_updates_include_zero_gradient_and_rollback` | the sign-consistent case gives `param == 855.2`, `param_lr == 4.8` |
| `compute_z(...)` | private `emg::tests::source_compute_z_selects_the_three_documented_regimes` | `z(mu - 1/60) == 0.471456263584609` |
| `emg_point(...)` | `emg_reference::source_scalar_model_goldens_cover_all_three_branches_and_units` | `emg_point(mu - 1/60) == 1992032.65711041` |
| `applyEstimatedParameters(...)` | `emg_reference::source_fixed_parameter_left_extension_has_exact_count_and_literal_endpoint` | 71 points with extension, first intensity `108845.941990663` |

Nothing is unaccounted. The three sections above five assertion macros -
`fitEMGPeakModel` for each container and `computeInitialMean` - are covered
against all seven of the source's literal traces, not only the cases the section
itself asserts.

Focused tests cover the source cutoff example, spectrum/chromatogram equality,
scalar model goldens, original sampling, extrapolation threshold crossings,
metadata and omitted arrays, explicit zero bounds, convergence and best-iterate
diagnostics, malformed input, numerical failure and shared resource limits.
Private unit tests independently check all four gradients by central
differences, training collection order, and iRprop+ sign/zero/rollback behavior.
Independent extracted reference tests cover all seven source traces, including
minute and second inputs. Fixture provenance records the original `f64`
decimals, the `f32` container conversions and source hashes.
