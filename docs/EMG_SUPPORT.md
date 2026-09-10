# Exponentially modified Gaussian fitting

`analysis::emg::EmgGradientDescent` ports the scientific operations of OpenMS
`MATH/MISC/EmgGradientDescent` from revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. It fits native spectra and
chromatograms, estimates parameters from `f64` slices, and evaluates supplied
parameters with optional extrapolation. It does not use C++ or an external
optimizer. `libm` 0.2.16 supplies only the complementary error function; the
remaining mathematical functions use Rust's standard library.

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

Focused tests cover the source cutoff example, spectrum/chromatogram equality,
scalar model goldens, original sampling, extrapolation threshold crossings,
metadata and omitted arrays, explicit zero bounds, convergence and best-iterate
diagnostics, malformed input, numerical failure and shared resource limits.
Private unit tests independently check all four gradients by central
differences, training collection order, and iRprop+ sign/zero/rollback behavior.
Independent extracted reference tests cover all seven source traces, including
minute and second inputs. Fixture provenance records the original `f64`
decimals, the `f32` container conversions and source hashes.
