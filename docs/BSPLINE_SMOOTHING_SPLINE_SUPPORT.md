# BSplineSmoothingSpline

Port of `src/openms/include/OpenMS/MATH/MISC/BSplineSmoothingSpline.h` and
`src/openms/source/MATH/MISC/BSplineSmoothingSpline.cpp` at openms4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Rust: `src/processing/spline/smoothing.rs`, re-exported as
`openms::processing::spline::BSplineSmoothingSpline`.

Tests: `tests/spline_math.rs` and the unit tests in
`src/processing/spline/smoothing.rs`.
Provenance: `tests/data/spline_math_provenance.json`.

## What the class actually does

The header says the class minimises `RSS(f) + lambda * integral(f''(x))^2 dx`,
or equivalently `integral(f''(x))^2` subject to `RSS(f) <= s`. It does not. It
exists to approximate `scipy.interpolate.UnivariateSpline` closely enough for
PyProphet, and its strategy is stated further down the same comment block: it
cannot modify the eol-bspline library, so it tries configurations and keeps the
one whose residual sum of squares lands nearest the budget.

Concretely, `s` is a *budget on the residual sum of squares*, and the algorithm
is:

1. Reject input that is not at least two strictly increasing points with
   matching ordinates.
2. Resolve the smoothing parameter: a negative `s` becomes scipy's default
   `m - sqrt(2m)`.
3. If there are fewer than four points, or the resolved `s` is zero or less, fit
   one `BSpline2d` with an automatic node count, report `n - 2` interior knots
   and stop — *without checking how close the fit is*.
4. Otherwise fit a cubic **polynomial** by normal equations. If its residual sum
   of squares is within ten percent of the budget, keep it and report zero
   interior knots.
5. Otherwise build a `BSpline2d` for each of the node counts `4, 6, 8, n/2,
   3n/4, n` (each raised to at least four, then sorted and deduplicated), keep
   those that fit, and choose the one whose residual sum of squares is closest to
   the budget, preferring fewer interior knots when two are within `0.001` of
   each other and both under budget.

So the "smoothing" comes from restricting the degrees of freedom, not from a
curvature penalty, and step 4 means a polynomial that misses the data can beat a
spline that matches it. Every one of those behaviours is reproduced.

## API mapping

Every public member of the header, in declaration order.

| C++ member | Rust | Notes |
|---|---|---|
| `BSplineSmoothingSpline(x, y, s = -1.0, k = 3)` | `BSplineSmoothingSpline::new(x, y)` for the defaults; `BSplineSmoothingSpline::with_smoothing(x, y, s, degree)` otherwise | Returns `Result`. Every path that sets `ok_ = false` in C++ — size mismatch, fewer than two points, non-increasing `x`, no usable candidate — is an `Err`. |
| `~BSplineSmoothingSpline()` | compiler-generated `Drop` | The C++ destructor is `= default`; it exists only so the `unique_ptr<BSpline2d>` can be destroyed against a complete type. Nothing to port. |
| `BSplineSmoothingSpline(const BSplineSmoothingSpline&) = delete` | `Clone` **is** implemented | The C++ deletes copying because it owns a `unique_ptr`. The Rust value owns its `BSpline2d` by value, so copying is cheap and safe; making it uncopyable would be a restriction with no cause. Deliberate divergence. |
| `operator=(const BSplineSmoothingSpline&) = delete` | covered by `Clone` | As above. |
| `BSplineSmoothingSpline(BSplineSmoothingSpline&&) = default` | Rust move semantics | Every Rust value moves by default. |
| `operator=(BSplineSmoothingSpline&&) = default` | Rust move semantics | As above. |
| `double eval(double x) const` | `BSplineSmoothingSpline::eval(&self, f64) -> Result<f64>` | The C++ returns `NaN` when `!ok_` or when a `BSPLINE` fit has a null pointer; neither state can exist here, because construction failure is an `Err`. |
| `bool ok() const` | **not ported as a method**; expressed as the `Result` from the constructors | A method that could only ever return `true` would be misleading. The class test's four `ok() == false` cases are asserted as `Err` instead. |
| `int num_interior_knots() const` | `BSplineSmoothingSpline::num_interior_knots(&self) -> i32` | Same bookkeeping, including that it counts grid nodes minus two rather than B-spline basis knots, and that the interpolating branch reports `n - 2` regardless of the grid it actually built. |
| `double rss() const` | `BSplineSmoothingSpline::rss(&self) -> f64` | Recomputed the same way, in the same order. |
| `double smoothing_param() const` | `BSplineSmoothingSpline::smoothing_param(&self) -> f64` | The resolved `s`, after a negative request has been replaced by `m - sqrt(2m)`. |
| `bool ok_`, `int num_interior_knots_`, `double rss_`, `double s_`, `int k_` (private) | private fields | `k_` is `[[maybe_unused]]` in the C++ and is surfaced as `degree()` here. |
| `std::vector<double> x_, y_` (private) | not stored | The C++ copies both inputs into members and never reads them again. Storing them would be dead weight; the fit itself keeps whatever it needs. |
| `enum FitType { POLYNOMIAL, BSPLINE }`, `fit_type_` (private) | private `enum Fit { Polynomial(Vec<f64>), Spline(Box<BSpline2d>) }` | Rust's enum carries the payload, so the tag and the two mutually exclusive members become one field. Surfaced read-only as `is_polynomial()`. |
| `std::vector<double> poly_coeffs_` (private) | the `Fit::Polynomial` payload | |
| `std::unique_ptr<BSpline2d> spline_` (private) | the `Fit::Spline` payload | |
| `void fit_smoothing_spline(...)` (private) | private `fit_smoothing_spline` | Ported including the candidate list and the comparator. |
| `bool try_polynomial_fit(...)` (private) | private `try_polynomial_fit` | Returns `Option<Vec<f64>>` rather than mutating members, which is what the C++ does anyway apart from the assignment at the end. |
| `double compute_rss(const BSpline2d*, ...)` (private) | private `spline_rss` | Returns `Err` rather than the C++ infinity for a null or not-ok spline, which cannot occur here. |
| `double compute_polynomial_rss(...)` (private) | private `polynomial_rss` | |
| `double eval_polynomial(double) const` (private) | private `eval_polynomial` | Ascending powers with a running product, in the source's order. |

Native additions:

| Rust | Why |
|---|---|
| `BSplineSmoothingSpline::degree()` | Reports the `k` the caller passed. The source stores it in a `[[maybe_unused]]` member and never reads it; the accessor exists so a caller can see what was asked for, and its documentation says plainly that it does not describe the curve. |
| `BSplineSmoothingSpline::is_polynomial()` | Read-only view of `fit_type_`, which decides whether extrapolation is cubic or decays to the fitted mean. |
| `BSplineSmoothingSpline::MAX_POINTS` | Bounded work; inherited from `BSpline2d::MAX_POINTS`. |

## Preserved source conventions

* **The scipy default is `m - sqrt(2m)` with `m` the point count**, computed as
  `m - std::sqrt(2.0 * m)`. For two points that is exactly zero, which then
  selects the interpolating branch — the probe's `smooth_n2` shows
  `smoothing_param = 0` for a negative request.
* **The ten-percent margin.** `if (rss <= s_target * 1.1)` accepts the
  polynomial. This is why `smooth_wiggle5_s2` reports a polynomial with residual
  `1.889` against a budget of `2.0` even though the spline branch would have
  produced `0.00375`.
* **The polynomial solve, exactly.** Normal equations `X^T X c = X^T y` built by
  accumulating ascending powers per point, Gaussian elimination with partial
  pivoting, a singularity threshold of `1e-10` on the pivot magnitude, and back
  substitution. Same loop order, same accumulation order.
* **The candidate node counts** `4, 6, 8, max(4, n/2), max(4, 3n/4), n` with the
  integer divisions as written, then `sort` and `unique`.
* **The candidate comparator, verbatim** — including that it is not a strict weak
  ordering. `std::sort` on at most six elements runs libstdc++'s insertion sort,
  which is stable, so the selection is well defined in practice; the port uses an
  explicit stable insertion sort with the same comparator rather than
  `sort_by`, so the agreement does not depend on a library's internal choice.
* **`k` is ignored.** The header's comment says the member is retained for object
  layout compatibility. The probe confirms it: `smooth_degrees_k1`, `k2` and `k3`
  produce identical residuals and identical values.
* **The interpolating branch does not check its fit.** For `n < 4` or `s <= 0` the
  constructor accepts whatever `BSpline2d` produces and reports `n - 2` interior
  knots, even though the grid it built has `2n` intervals. Reproduced.
* **"Interpolation" does not interpolate.** With `s = 0` the underlying fit is
  still least squares; the probe's `smooth_linear5_zero` misses the first point
  by `7.8e-11` and `smooth_n2` misses it by `0.041`. The upstream class test
  allows `0.05` for exactly this reason.

## Native differences

| Difference | Source behaviour | Port behaviour | Why |
|---|---|---|---|
| Failed fit | `ok()` false; `eval` returns `NaN`; `num_interior_knots`, `rss` and `smoothing_param` return whatever the partially-run constructor left. | `Err(Error::InvalidValue)`, or `Error::UnsortedData` for non-increasing `x`. | A value that only reports `NaN` is not a value. The four `ok() == false` cases in the class test map one-to-one onto the four `Err` cases. |
| `ok()` | Public method. | Not present; the constructor's `Result` carries it. | See above. |
| Non-finite input or `s` | Not checked. A `NaN` `s` fails `s < 0.0`, so it is used as the budget and every comparison against it is false. | `Error::InvalidValue`. | No silent `NaN` propagation. |
| Copying | Deleted. | `Clone` implemented. | The deletion exists only because of the `unique_ptr`. |
| Error logging | `OPENMS_LOG_ERROR` on each failure path. | The `Err` message names the condition. | The crate does not log; the message carries the same information to the caller instead of to a stream. |
| Point count | Unbounded. | `MAX_POINTS`, checked before anything is allocated. | Bounded work. |
| Candidate sort | `std::sort` with a comparator that is not a strict weak ordering. | Explicit stable insertion sort with the same comparator. | Same result for the at-most-six candidates libstdc++ insertion-sorts, without depending on that implementation detail. |

No OpenMP in `BSplineSmoothingSpline.cpp`; nothing to record.

## Checked boundaries and evidence

Evidence tier 2 (executed probe), bit-for-bit, against 27 fits of the
unmodified C++ covering every branch: the polynomial branch, the interpolating
branch at `n = 2`, `3`, `5` and `6`, the node-count search, and the four failure
paths. Flags and hashes are in `tests/data/spline_math_provenance.json`; the
contraction caveat in `docs/BSPLINE2D_SUPPORT.md` applies here too, since this
class is built on `BSpline2d`.

### Class-test sections

`BSplineSmoothingSpline_test.cpp` has 16 `START_SECTION`s and 53 assertion
macros. All 16 are mapped. Unless stated otherwise the Rust test is
`smoothing_spline_reproduces_every_probe_case` in `tests/spline_math.rs`, which
compares `num_interior_knots`, `rss`, `smoothing_param` and every sampled value
for the named probe case.

| Section | Probe case / Rust test | One value it reproduces |
|---|---|---|
| `BSplineSmoothingSpline(const std::vector<double>&, const std::vector<double>&, double, int)` | `smooth_linear5_auto`, `smooth_linear5_zero`, `smooth_linear5_s2`; unit test `the_scipy_default_budget_selects_a_cubic_polynomial` | `smoothing_param() == 1.8377223398316205` for five points, which is `5 - sqrt(10)` |
| `~BSplineSmoothingSpline()` | no Rust construct | The C++ destructor is `= default`. The Rust value owns its fit and has no `Drop` impl; every test drops many of them. Mapped as a native equivalent with no behaviour to assert. |
| `double eval(double) const` | `smooth_square5_auto` | `eval(1.5) == 2.25` on `y = x^2`, and the section's extrapolation checks: `eval(-0.5) == 0.25` and `eval(4.5) == 20.25` are finite |
| `bool ok() const` | `smooth_ramp5_auto` for the valid case; `smooth_size_mismatch`, `smooth_single_point`, `smooth_unsorted` for the three invalid ones | The three invalid constructions return `Err`, matching the probe's `ok = 0`; the valid one gives `eval(2.0) == 2.9999999999999947` |
| `int num_interior_knots() const` | `smooth_ramp6_zero`, `smooth_ramp6_s10` | The interpolating branch reports `4` for six points, the smoothing branch `0` |
| `double rss() const` | `smooth_ramp5_auto` | `rss() == 1.5399057907979914e-27`, finite and non-negative as the section asserts |
| `double smoothing_param() const` | `smooth_ramp6_auto`, `smooth_ramp6_s5`, `smooth_ramp6_zero` | `smoothing_param() == 2.5358983848622456` for six points, inside the section's `(2, 3)` window, and `5.0` and `0.0` are returned verbatim |
| `Smoothing behavior test` | `smooth_noisy20_s20` | Twenty sine points with a deterministic sawtooth perturbation and a budget of 20: `rss() == 0.52189567174737794`, `eval(1.5) == 1.0118796941908261` |
| `Interpolation vs smoothing test` | `smooth_wiggle5_zero`, `smooth_wiggle5_s2`; unit test `a_polynomial_within_the_ten_percent_margin_wins_over_a_closer_spline` | Interpolation gives `eval(1.0) == 2.9685980812295538`, within the section's `0.05` of the data point `3.0`; smoothing gives `eval(0.0) == 1.1642857142857219` |
| `Linear data test` | `smooth_line6_auto` | `eval(3.0) == 7.0000000000000124` against the exact `2*3+1`, and `rss() == 3.3372760595374904e-27`, well under the section's `0.1` |
| `Quadratic data test` | `smooth_quad6_s1` | `eval(-1.0) == 1.0000000000000007`, within the section's `1.0` of `1.0` |
| `Edge case: identical x values should fail` | `smooth_duplicate_x` | Probe `ok = 0`; the port returns `Error::UnsortedData` |
| `Edge case: small dataset` | `smooth_n2`; unit test `small_datasets_fall_back_to_the_interpolating_branch` | `eval(0.5) == 0.5`, and `eval(0.0) == 0.040851301587258748` — the fit does not pass through the data |
| `Edge case: n=3 dataset` | `smooth_n3` | `num_interior_knots() == 1` and `eval(2.0) == 2.9999999998164566` |
| `Different spline degrees` | `smooth_degrees_k1`, `k2`, `k3`; unit test `the_degree_argument_changes_nothing` | All three give `rss() == 0.96230158730158666` and `eval(2.0) == 2.0793650793650973` |
| `Consistency test: repeated construction` | unit test `repeated_construction_is_deterministic` | Two `with_smoothing(x, y, 2.0, 3)` values agree exactly on `rss()`, `num_interior_knots()` and every sampled value |

Unaccounted sections: none.

### Beyond the class test

The class test never reaches the node-count search: every one of its datasets is
either small enough for the interpolating branch or close enough to a cubic for
the polynomial branch. The probe adds an alternating `0, 1, 0, 1, ...` sequence
of twelve points, for which no cubic gets inside a budget of `0.5`; the search
then runs, builds five candidates and selects the densest grid, reporting
`num_interior_knots() == 10` and `rss() == 0.0018445489048118143`. The unit test
`a_budget_no_polynomial_can_meet_runs_the_node_search` covers it, together with
the same data at a budget of `3.0` where the polynomial wins again.

### Independently derived checks

* `smoothing_param()` for a negative request equals `m - sqrt(2m)` for the point
  count, checked as arithmetic rather than transcribed.
* Two points make the scipy default exactly zero, which is why `n = 2` always
  takes the interpolating branch regardless of the `s` requested.
* Far outside the node domain a `BSPLINE` fit returns the fitted mean: the
  alternating twelve-point case gives `eval(13.0) == 0.5`, which is the mean of
  six zeros and six ones.
* Repeated construction from the same input is bit-identical, which is the
  section's own claim and is not something the probe can establish.

### Boundaries the port checks

* At least two points, matching lengths, at most `MAX_POINTS`.
* All values and the smoothing parameter finite.
* Strictly increasing abscissae (`Error::UnsortedData`).
* At least one usable candidate in the node search.
* Every `BSpline2d` the fit builds carries that type's own ceilings and
  finiteness checks.
