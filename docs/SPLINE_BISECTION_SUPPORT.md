# SplineBisection

Port of `src/openms/include/OpenMS/MATH/MISC/SplineBisection.h` at openms4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. Header-only: the header *is* the
implementation, a single 34-line function template.

Rust: `src/processing/spline/bisection.rs`, re-exported as
`openms::processing::spline::{spline_bisection, SplineFunction}`.

Tests: `tests/spline_math.rs` and the unit tests in
`src/processing/spline/bisection.rs`.
Provenance: `tests/data/spline_math_provenance.json`.

## API mapping

The header has one public entity.

| C++ member | Rust | Notes |
|---|---|---|
| `template <class T> void spline_bisection(const T& peak_spline, double const left_neighbor_mz, double const right_neighbor_mz, double& max_peak_mz, double& max_peak_int, double const threshold = 1e-6)` | `pub fn spline_bisection<T: SplineFunction + ?Sized>(spline: &T, left_neighbor: f64, right_neighbor: f64, threshold: f64) -> Result<(f64, f64)>` | The two out-parameters become the return tuple `(position, value)`. Rust has no default arguments, so the source's `1e-6` is exported as `DEFAULT_BISECTION_THRESHOLD` and passed explicitly. |
| the implicit requirement `T::eval(double) const` and `T::derivative(double) const` | `trait SplineFunction { fn eval(&self, f64) -> Result<f64>; fn first_derivative(&self, f64) -> Result<f64>; }` | The template's duck typing becomes an explicit bound. `first_derivative` is not called `derivative` because `CubicSpline2d::derivative` takes a derivative order, carrying the C++ `derivatives` overload. Implemented for both types the header's comment names: `CubicSpline2d` and `BSpline2d`. |

Native additions:

| Rust | Why |
|---|---|
| `MAX_BISECTION_STEPS` | The source's loop can run forever; see below. |
| `DEFAULT_BISECTION_THRESHOLD` | The source's default argument, made nameable. |

Related but distinct: `CubicSpline2d::peak_maximum(left, right, tolerance)` is a
native function that predates this module and is used by
`processing::peak_picking`. It differs deliberately — it validates that both
bracket ends are inside the spline's domain, it also stops when the midpoint
stops moving, and it fails after 128 halvings — and is *not* a port of this
header. Use `spline_bisection` when source behaviour is what is wanted.

## Preserved source conventions

The header is short, and almost all of it is the termination condition, so the
port is a transcription with one addition. Item by item:

1. **`do` / `while`, not `while`.** The body always runs at least once. With
   `right <= left` the width test is false from the start, so exactly one
   midpoint is evaluated and the function returns the midpoint of the *original*
   bracket. The probe's `reversed_bracket` case passes `(501, 499)` and gets
   `500.0` back.
2. **`lefthand_sign` is seeded `true` and never updated.** The source computes
   `lefthand_sign ^ midpoint_sign`, which therefore reduces to
   `!midpoint_sign`: move the right end down when the derivative is negative,
   move the left end up otherwise. The port keeps the reduced form with a comment
   naming what it reduces from, because writing the exclusive-or against a
   constant would be theatre.
3. **`midpoint_sign = (d < 0.0) ? false : true`.** A derivative of exactly `-0.0`
   compares `>= 0.0` and counts as positive. The Rust spelling
   `midpoint_deriv_val >= 0.0` has the same behaviour on `-0.0`.
4. **The early exit is spelled `!(fabs(d) > eps)` with `eps = DBL_EPSILON`.**
   That is not the same as `fabs(d) <= eps`: it is also true when the derivative
   is `NaN`. The port writes both halves out —
   `magnitude <= eps || magnitude.is_nan()` — so the `NaN` case survives a
   future edit, and so clippy does not rewrite the negation.
5. **The reported position after an early exit is not the midpoint that
   triggered it.** The loop breaks *before* narrowing, then reports
   `(lefthand + righthand) / 2`. When the break happens on the first pass those
   are still the caller's bracket ends.
6. **The threshold is compared against the bracket width, not the step.**
   `while (righthand - lefthand > threshold)`.
7. **The apex may lie outside the bracket.** The derivative then never changes
   sign and the search converges onto the bracket end nearest the true apex. The
   class test pins the right-hand case; the probe adds the left-hand one and the
   case where the apex sits exactly on the left edge.
8. **The value is `eval` at the reported position**, computed after the loop.

## Native differences

| Difference | Source behaviour | Port behaviour | Why |
|---|---|---|---|
| Iteration bound | None. The loop halves until the width test fails. When the bracket ends are adjacent `f64` values and their exact midpoint ties to the end the derivative's sign selects, the bracket stops shrinking and the loop never terminates. A `threshold` of zero or a negative one makes that certain. | `MAX_BISECTION_STEPS = 4096`, then `Error::InvalidValue`. | No hangs on untrusted input. A bisection across the whole `f64` range needs about 2 100 halvings, so a search the C++ would finish never reaches the ceiling. The unit test `a_threshold_below_the_float_spacing_hits_the_ceiling_instead_of_hanging` constructs the non-terminating case exactly: bracket `(1.0, 1.0 + 2^-52)`, whose midpoint `1 + 2^-53` is a tie that rounds to the even mantissa and so back to the left end, with an apex far to the right to keep the derivative positive. |
| Threshold validation | None; zero, negative and `NaN` are accepted and all defeat the loop. | Must be finite and strictly positive. | Same reason. |
| Bracket validation | None. | Both ends must be finite. | `(NaN + x) / 2` is `NaN`, and the derivative of a `NaN` is `NaN`, which the early exit then swallows — returning a `NaN` position with a straight face. |
| Errors from the spline | `eval` and `derivative` may throw; the exception propagates out of the template. | `?` propagates the `Err`. | Same shape. A `CubicSpline2d` bracket outside the knot range is the common case, and `spline_bisection_over_a_cubic_spline_finds_the_probe_apex` asserts it. |
| Out-parameters | Two `double&`. | Returned tuple `(position, value)`. | Rust convention; the values and their order are unchanged. |

No OpenMP; the header has no pragma.

## Checked boundaries and evidence

Evidence tier 2 (executed probe), bit-for-bit. The probe compiles the header
unmodified and drives it with the class test's own `ParabolaSpline` model, then
with a real `CubicSpline2d` over the peak from `CubicSpline2d_test.cpp`.

### Class-test sections

`SplineBisection_test.cpp` has 1 `START_SECTION` and 7 assertion macros. It is
mapped.

| Section | Rust test | One value it reproduces |
|---|---|---|
| `template <class T> void spline_bisection(const T&, double, double, double&, double&, double)` | `spline_bisection_reproduces_every_probe_bracket`, and the unit test `reproduces_the_probe_brackets` | The section's four cases in order: the centred apex gives `(500.0, 1000.0)`; the off-centre one `(500.29999971389771, 749.99999999999989)`; the tighter threshold `(500.12345599988475, 100.0)`; and the apex outside the bracket converges to `500.99999952316284`, which is the section's `TEST_REAL_SIMILAR(mz, 501.0)` at its `1e-4` tolerance. All compared with exact equality here. |

Unaccounted sections: none.

### Beyond the class test

The probe adds four cases the section does not cover, each pinning a branch the
prose above describes:

* apex outside the bracket on the *left*, converging to `499.00000047683716`;
* apex exactly on the left edge, giving `(499.00000047683716, 9.9999999999997726)`;
* a reversed bracket `(501, 499)`, where the `do`/`while` runs once and returns
  `500.0`;
* a coarse threshold of `0.1`, which stops at `500.28125` instead of `500.3`.

It also runs the search over a `CubicSpline2d` — the shape the header is
documented for — bracketing the upstream peak between `486.795` and `486.800`
and finding `486.79738311767579` at the default threshold and
`486.79738342076541` at `1e-9`.

### Independently derived checks

* The parabola model is exact, so the apex of `-a(x - p)^2 + h` is `p` with value
  `h`; the centred and off-centre cases converge to those within the threshold
  they were given, which is a property of bisection rather than a transcribed
  number.
* A tighter threshold gives a strictly closer answer: `1e-9` reaches
  `500.12345599988475` where `1e-6` would not.
* Driving a `BSpline2d` far outside its fitted domain, where the derivative is
  exactly zero, exits on the first pass and returns the bracket midpoint —
  `spline_bisection_also_drives_a_bspline` asserts `-50.0` for the bracket
  `(-60, -40)`, which is the `do`/`while`'s first-pass behaviour.

### Boundaries the port checks

* Both bracket ends finite.
* `threshold` finite and strictly positive.
* At most `MAX_BISECTION_STEPS` halvings.
* Errors from the driven spline propagate unchanged; nothing is swallowed.
