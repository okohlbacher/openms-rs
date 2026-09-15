# CubicSpline2d

Port of `src/openms/include/OpenMS/MATH/MISC/CubicSpline2d.h` and
`src/openms/source/MATH/MISC/CubicSpline2d.cpp` at openms4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Rust: `src/processing/spline/cubic.rs`, re-exported as
`openms::processing::spline::CubicSpline2d` and, for the callers that predate
this module, as `openms::processing::peak_picking::CubicSpline2d`.

Tests: `tests/spline_math.rs` and the unit tests in `src/processing/spline/cubic.rs`.
Fixtures: `tests/data/spline_math_cpp_probe.tsv` (the C++ oracle) and
`tests/data/cubic_spline_picker_supports.tsv` (real picker inputs, port data
rather than oracle data).
Provenance: `tests/data/spline_math_provenance.json`.

## API mapping

Every public member of the header, in declaration order.

| C++ member | Rust | Notes |
|---|---|---|
| `CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y)` | `CubicSpline2d::new(&[f64], &[f64]) -> Result<Self>` | Also `with_max_points(x, y, max_points)` for an explicit knot ceiling. The three `Exception::IllegalArgument` conditions become `Error::InvalidValue`; two further conditions the C++ does not check are added (below). |
| `CubicSpline2d(const std::map<double, double>& m)` | `CubicSpline2d::from_pairs(&[(f64, f64)]) -> Result<Self>` | The map's sorting and key deduplication are reproduced explicitly: pairs are sorted by abscissa and, among equal abscissae, the first ordinate wins, as `std::map::insert` does. The "fewer than two entries" check therefore applies to the *deduplicated* count, as in the C++. |
| `double eval(double x) const` | `CubicSpline2d::eval(f64) -> Result<f64>` | Same segment selection, same Horner form, same evaluation order. |
| `double derivative(double x) const` | `CubicSpline2d::derivative(f64, 1)` | The C++ body is `return derivatives(x, 1);`, so the two overloads collapse into one Rust method with an order argument. |
| `double derivatives(double x, unsigned order) const` | `CubicSpline2d::derivative(f64, u8) -> Result<f64>` | `order` is `u8`; 0 and anything above 3 are rejected, as `order < 1 \|\| order > 3` is in the C++. The range check runs before the order check, as in the C++. |
| `void init_(...)` (private) | the body of `with_max_points` | Not public API; listed because it is where the recurrence lives. |
| `a_`, `b_`, `c_`, `d_`, `x_` (private) | private fields of the same name | Not public API. `domain()` and `segment_count()` are native accessors added because the peak picker needs the knot range. |

Native additions, all documented at the item:

| Rust | Why |
|---|---|
| `CubicSpline2d::MAX_POINTS` and `with_max_points` | The source allocates five vectors proportional to the knot count with no ceiling. |
| `CubicSpline2d::domain()` | The knot range; needed by every caller that must decide whether a query is legal before making it. |
| `CubicSpline2d::segment_count()` | Number of cubic segments; used by the map-constructor test to prove deduplication happened. |
| `CubicSpline2d::peak_maximum(left, right, tolerance)` | Native predecessor of `spline_bisection`, retained because `processing::peak_picking` depends on its extra guarantees (both bracket ends validated, stops when the midpoint stops moving, fails after 128 halvings). `docs/SPLINE_BISECTION_SUPPORT.md` compares the two. |
| `impl SplineFunction for CubicSpline2d` | Makes the type usable with `spline_bisection`, which is what the C++ template's duck typing achieves. |
| `CubicSpline2dFitter` | Reusable construction storage for a caller that fits many splines in a loop. Same recurrence, same coefficients, one set of buffers instead of eight allocations per spline. See *Reusable construction storage* below. |

## Reusable construction storage

`with_max_points` allocates eight vectors per spline — `h`, `mu` and `z` for the
tridiagonal sweep, `b`, `c` and `d` for the coefficients, a copy of the knots and
a copy of the leading ordinates — and frees all eight when the spline is dropped.
For a caller that wants one spline that is invisible. `processing::peak_picking`
fits one spline per candidate centroid over a handful of knots and drops it
before the next peak, so for it those eight allocations are the construction's
largest single item: the profiling lane measured `alloc::alloc` at 271 144 136 of
the 590 918 793 instructions `with_max_points` costs, 46 %.

`CubicSpline2dFitter` owns the eight buffers instead. `fit` clears and refills
them and hands the spline back by reference, so a loop pays for allocation only
while the buffers grow to the largest support it has seen. The arithmetic is not
touched: both entry points call one private `fit_into`, which is the transcribed
recurrence, so there is no second copy of it that could drift.

### Measured

Benchmark input `UK222.mzML` (2.3 GB, 40 856 spectra, 197 765 338 raw points), the
INI the C++ `PeakPickerHiRes` wrote, one thread, pinned, on `ibminode06`. The run
constructs **13 856 120** splines, so eight allocations per spline is
**110 848 960** malloc/free pairs.

Bit-identity was checked at that scale, not sampled: every one of the 13 856 120
splines was fitted by `with_max_points` before the change, by `with_max_points`
after it and by `CubicSpline2dFitter`, and all knots and all four coefficient
vectors of each were written out and hashed. The three streams are the same
4 580 021 720 bytes, sha256
`53d8531925a4f83e52c409307769830ab8f4edb5ecbebf834117833b664c91c8`. The end-to-end
tool output is unchanged as well:
`bb13eecfe092a272b08ddc71feec3780c7bc876e8847e8a45b173cda9d2dad52`, 535 613 726
bytes, with and without the fitter. `tests/data/cubic_spline_picker_supports.tsv`
keeps a stratified sample of those supports so the unit test
`a_reused_fitter_replays_real_picker_supports_coefficient_for_coefficient`
re-checks the property on real data in every run of the suite.

Instruction counts are callgrind simulations, which do not depend on node load.
The harness replays the recorded supports and does per record what the picker's
inner loop does around the spline: construct, then bisect for the maximum. The
numbers below are for the 214 780 supports of the profiling lane's 682-spectrum
slice, whole program, so they include a fixed parse and bisection cost that is
the same in all three columns.

| | instructions | vs. before |
|---|---:|---:|
| `with_max_points`, before this change | 1 331 496 468 | — |
| `with_max_points`, after | 1 217 243 728 | −8.6 % |
| `CubicSpline2dFitter::fit` | 916 455 502 | −31.2 % |

The one-shot path got faster too, which was not the point but is worth
recording: sizing the eight vectors with `Vec::with_capacity` and filling them
beats `vec![0.0; n]` plus `collect` plus `to_vec`, because `calloc` on a small
block memsets anyway and costs more to reach.

In wall time the fitter is worth about 1.3 s of the tool's 33 s, measured with
`PeakPickerHiRes` on the input above, one thread, pinned to an idle core pair
on `ibminode06`, eleven interleaved A/B repetitions: 33.19 s median before,
32.01 s after, paired difference 1.27 s median (user time 1.31 s). A sampled
phase trace puts all of it in the pick phase — 17.24 s to 15.11 s — with load
and write unchanged. This branch on its own, with the picker still calling
`with_max_points`, is 0.06 s slower and 0.09 s less CPU than before, i.e. the
same within run-to-run spread; the saving needs the picker's call site to move.

The picker's own call site is one line: `CubicSpline2d::with_max_points(&support.xs,
&support.ys, self.max_points)?` becomes `fitter.fit_with_max_points(&support.xs,
&support.ys, self.max_points)?` with `let mut fitter = CubicSpline2dFitter::new();`
hoisted next to `let mut support = Support::default();`. Nothing downstream
changes: `spline_bisection(&spline, ..)` still compiles because
`SplineFunction` is implemented for `&T` (see `docs/SPLINE_BISECTION_SUPPORT.md`),
and `spline.eval(..)` and `half_height(&spline, ..)` reach the value through the
ordinary deref coercion.

## Preserved source conventions

* **The recurrence, term by term and in order.** `init_` is transcribed with its
  original groupings: `l = 2 * (x[i+1] - x[i-1]) - h[i-1] * mu[i-1]`,
  `z[i] = (3 * (y[i+1]*h[i-1] - y[i]*(x[i+1]-x[i-1]) + y[i-1]*h[i]) / (h[i-1]*h[i]) - h[i-1]*z[i-1]) / l`,
  and the backward sweep `c[j] = z[j] - mu[j]*c[j+1]`,
  `b[j] = (y[j+1]-y[j])/h[j] - h[j]*(c[j+1] + 2*c[j])/3`,
  `d[j] = (c[j+1]-c[j])/(3*h[j])`. No term is factored, reassociated or
  precomputed, because floating-point addition is not associative and the probe
  compares bit-for-bit.
* **The natural boundary condition.** Confirmed from the source, not assumed:
  `CubicSpline2d.cpp:152` sets `c_.back() = 0`, and the forward sweep never
  writes `mu[0]` or `z[0]`, so `c_[0] = z[0] - mu[0]*c_[1] = 0`. The second
  derivative is `2*c_[i] + 6*d_[i]*t`.

  The two ends differ, and the difference is measurable rather than cosmetic. At
  the first knot `t` is zero, so the value is `2*c_[0]`, which is **exactly**
  `0.0` for every input. At the last knot the evaluation falls on the last
  segment at its right end, `2*c_[n-1] + 6*d_[n-1]*h` with
  `d_[n-1] = (c_[n] - c_[n-1]) / (3*h)`. Algebraically that is `2*c_[n] = 0`,
  but three roundings stand between the two, so it is zero only up to rounding.
  The probe shows both outcomes: `cubic_sine d2_last` is exactly `0` on the
  class test's uniform grid, while `cubic_upstream d2@486.811` is
  `-3.814697265625e-06` against interior second derivatives of order `1e11` —
  about `3e-17` relative. The C++ behaves identically, and its own class test
  uses `TEST_REAL_SIMILAR(sp5.derivatives(x[n], 2), 0)` rather than an equality.
  The Rust test `the_natural_boundary_condition_holds_at_both_ends` asserts
  exact `0.0` at the first knot of both fixtures and the probe's non-zero value
  at the last knot of the peak, so the distinction cannot silently rot.
* **Segment selection at a knot.** `lower_bound` followed by
  `if (x_[i] > x || x_.back() == x) --i` picks the segment *starting* at an
  interior knot, so `eval(x_i)` returns `y_i` with `t = 0`; at the last knot it
  steps back one segment and evaluates at its right end. `partition_point`
  produces the same index for strictly increasing abscissae.
* **Out-of-domain queries are errors.** `x < x_.front() || x > x_.back()` throws
  in the C++ and returns `Error::InvalidValue` here. This spline does not
  extrapolate — unlike `BSpline2d`, which does; see
  `docs/BSPLINE2D_SUPPORT.md`.
* **Derivative order 3 is discontinuous at the knots.** Kept, because it follows
  from the piecewise-cubic form: `derivative(x, 3)` is `6*d_[i]`, constant on
  each segment. Orders 1 and 2 are continuous, which the Rust test checks across
  every interior knot.

## Native differences

| Difference | Source behaviour | Port behaviour | Why |
|---|---|---|---|
| Repeated abscissae | `adjacent_find(x.begin(), x.end(), std::greater<double>())` only rejects a *decrease*, so `x = {0, 1, 1, 2}` is accepted; `h[i] = 0` then divides in `mu[i] = h[i]/l` and `d[j] = (c[j+1]-c[j])/(3*h[j])`, and the spline evaluates to `NaN`. The probe records `cubic_duplicate_x threw=0, eval_is_nan=1`. | `Error::InvalidValue` at construction. | A `NaN` spline is reported by nothing and propagates into the caller's data. Rejecting costs one pass over `x`. `from_pairs` is the exception and keeps the source's meaning: there a repeat is a single knot, because `std::map` has already collapsed it. |
| Non-finite input | Not checked. `NaN` defeats `adjacent_find`, so a `NaN` abscissa passes and poisons the coefficients. | `Error::InvalidValue`. | Same reasoning. |
| Non-finite intermediates and results | Not checked; an overflowing knot spacing gives `inf` coefficients. | Every recurrence intermediate and every evaluation goes through a finiteness check and yields `Error::InvalidValue`. | No silent `NaN`/`inf` propagation. |
| Non-finite query position | `x < front()` and `x > back()` are both false for `NaN`, so the range check passes and `lower_bound` returns `end()`, which the C++ then indexes. | `Error::InvalidValue`. | The C++ path reads one past the end of `x_`. |
| Knot count | Unbounded. | `MAX_POINTS = 1 000 000` by default, overridable per call. | Bounded work; checked before anything is allocated. |
| `derivative` / `derivatives` | Two overloads. | One method with an order argument. | Rust has no overloading; the C++ `derivative` body is literally `derivatives(x, 1)`. |

No OpenMP: `CubicSpline2d.cpp` carries no `#pragma omp`, so there is no
parallelism gap to record.

## Checked boundaries and evidence

Evidence tier 2 (executed probe). `tests/data/spline_math_cpp_probe.tsv` is the
output of a driver that compiles `CubicSpline2d.cpp` unmodified against shims
for the OpenMS build-configuration headers and prints values at `%.17g`; the
driver, its flags and its sha256 are in `tests/data/spline_math_provenance.json`.
Every comparison below is bit-for-bit equality, not a tolerance.

The class test's own literals (tier 3) agree with the probe to every digit they
print: `eval(486.785) = 35173.1841778984`, `eval(486.794) = 2271426.93316241`,
`derivatives(486.785, 1) = 39270152.2996247`,
`derivatives(486.785, 2) = 12290904368.2736`,
`derivatives(486.794, 1) = 594825947.154264`,
`derivatives(486.794, 2) = 7415503644.8958`.

### Class-test sections

`CubicSpline2d_test.cpp` has 4 `START_SECTION`s and 30 assertion macros. All 4
are mapped.

| Section | Rust test | One value it reproduces |
|---|---|---|
| `CubicSpline2d(const std::vector<double>&, const std::vector<double>&)` | `cubic_spline_reproduces_every_probe_row`, and the unit test `matches_the_cpp_probe_bit_for_bit` | `eval(486.785) == 35173.184177898438` |
| `CubicSpline2d(const std::map<double, double>&)` | `cubic_spline_reproduces_every_probe_row` (case `cubic_upstream_map`) and the unit test `the_map_constructor_sorts_and_keeps_the_first_ordinate_per_abscissa` | `eval(486.790) == 620386.5` from a reversed, repeat-bearing pair list |
| `double eval(double x)` | `cubic_spline_reproduces_every_probe_row` | `eval(486.794) == 2271426.9331624084` |
| `double derivatives(double x, unsigned order)` | `cubic_spline_reproduces_every_probe_row`, plus `the_natural_boundary_condition_holds_at_both_ends` and `first_and_second_derivatives_are_continuous_across_every_knot` | `derivative(486.785, 2) == 12290904368.273579`, and `derivative(x_first, 2) == 0.0` exactly |

Unaccounted sections: none.

### Independently derived checks

These do not come from the probe and would survive it being wrong:

* `derivative(x_first, 2)` is compared with exact `0.0`, which the recurrence
  forces bit for bit on every input: `c_[0]` is never written and the offset into
  the first segment is zero. The same is *not* claimed at the last knot, where
  the algebraic zero passes through three roundings; that end is pinned against
  the probe instead.
* The first and second derivatives are compared just left and right of every
  interior knot; a cubic spline is `C^2`, so a discontinuity there would mean the
  tridiagonal solve is wrong even if every probe value matched.
* `eval` at each input knot returns that knot's ordinate.
* `from_pairs` on a reversed pair list with an added repeat produces a spline
  with ten segments and the same values as the ordered construction, which is
  the defining property of the map constructor.

### Boundaries the port checks

* Matching slice lengths, at least two knots, at most `MAX_POINTS`. This is the
  only check that runs before `with_max_points` sizes its eight buffers, and it
  is the one that bounds them; the two below run after. An input that is
  rejected for a non-finite or a non-increasing abscissa therefore allocates and
  frees at most `max_points` knots' worth first, where the revision before the
  fitter allocated nothing. The work stays bounded by the same ceiling.
* All abscissae and ordinates finite, and the abscissae strictly increasing.
* Every recurrence intermediate finite, so an error leaves no half-built spline:
  `with_max_points` builds the coefficients into a value that is returned only
  after the last check passes and is dropped otherwise, and
  `CubicSpline2dFitter::fit` returns `Err` without handing out a reference to
  the spline it was filling in, so a half-written fitter is not observable — the
  next fit rewrites every buffer.
* Queries inside the closed knot range and finite; derivative order in 1..=3.
* `peak_maximum` validates both bracket ends and its tolerance, and fails after
  128 halvings rather than looping.
