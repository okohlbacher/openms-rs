# BSpline2d

Port of `src/openms/include/OpenMS/MATH/MISC/BSpline2d.h` and
`src/openms/source/MATH/MISC/BSpline2d.cpp` at openms4-core
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`.

Rust: `src/processing/spline/b_spline.rs`, re-exported as
`openms::processing::spline::BSpline2d`.

Tests: `tests/spline_math.rs` and the unit tests in
`src/processing/spline/b_spline.rs`.
Provenance: `tests/data/spline_math_provenance.json`.

## This is a reimplementation, not a wrapper — read this first

`BSpline2d.cpp` is 60 lines. It is a PIMPL wrapper: the constructor news an
`eol_bspline::BSpline<double>`, and `solve`, `eval`, `derivative` and `ok`
forward to it. The algorithm lives in the vendored `eol-bspline` library at
`src/openms/extern/eol-bspline/BSpline/` (UCAR, BSD-3-Clause), about 1 100 lines
of template C++ implementing Ooyama's cubic B-spline (Monthly Weather Review
115, October 1987).

This crate takes no third-party dependency, so **the library was reimplemented
in Rust from that pinned vendored source**, not wrapped and not approximated.
Concretely, the following are ported member for member into
`src/processing/spline/b_spline.rs`:

| eol-bspline | Rust |
|---|---|
| `BSplineBase<T>::Setup` | `setup` |
| `BSplineBase<T>::Ratiod` | `ratios` |
| `BSplineBase<T>::Alpha` | the tail of `setup` |
| `BSplineBase<T>::Beta` | `Domain::beta` |
| `BSplineBase<T>::Basis` | `Domain::basis` |
| `BSplineBase<T>::DBasis` | `Domain::dbasis` |
| `BSplineBase<T>::qDelta` | `Domain::q_delta` |
| `BSplineBase<T>::calculateQ` | `Domain::calculate_q` |
| `BSplineBase<T>::addP` | `Domain::add_p` |
| `BandedMatrix<T>` | `Band` |
| `LU_factor_banded` | `lu_factor_banded` |
| `LU_solve_banded` | `lu_solve_banded` |
| `BSpline<T>::solve` | `BSpline2d::solve` |
| `BSpline<T>::evaluate` | `BSpline2d::eval` |
| `BSpline<T>::slope` | `BSpline2d::derivative` |
| `BSpline<T>::coefficient` | `BSpline2d::coefficient` |

What that means for the evidence: a reimplementation can be wrong in ways a
wrapper cannot, so the whole surface is checked against an **executed probe** of
the unmodified C++ (tier 2), including the internal state the OpenMS wrapper
hides. For the 202-point class-test fixture the port reproduces all **405** basis
coefficients, the node count, the node domain, the penalty weight and every
sampled value and slope **bit-for-bit**. That is the claim this document rests
on; it is not a claim of API-level similarity.

## API mapping

Every public member of `BSpline2d.h`, in declaration order.

| C++ member | Rust | Notes |
|---|---|---|
| `enum BoundaryCondition { BC_ZERO_ENDPOINTS = 0, BC_ZERO_FIRST = 1, BC_ZERO_SECOND = 2 }` | `enum BoundaryCondition { ZeroEndpoints = 0, ZeroFirst = 1, ZeroSecond = 2 }` | Discriminants kept. The header's "don't change these, they are passed through to the eol-bspline implementation" no longer describes an external boundary, but the values still index the beta table and a caller may have stored one. `ZeroSecond` is `Default`, matching the C++ default argument. |
| `BSpline2d(x, y, wavelength = 0, boundary_condition = BC_ZERO_SECOND, num_nodes = 0)` | `BSpline2d::new(x, y)` for the defaults; `BSpline2d::with_options(x, y, wavelength, boundary_condition, num_nodes)` for the rest | Rust has no default arguments. Returns `Result`: every state in which the C++ would leave `ok()` false is an `Err`. |
| `virtual ~BSpline2d()` | compiler-generated `Drop` | The C++ destructor exists only to `delete spline_`. The Rust value owns `Vec<f64>`s and needs no manual drop; there is no `impl Drop`. |
| `bool solve(const std::vector<double>& y)` | `BSpline2d::solve(&mut self, &[f64]) -> Result<()>` | Refits new ordinates over the same domain. The C++ `OPENMS_PRECONDITION(nX() == y.size())` is inert in release builds and then reads past the end of `y`; here a length mismatch is `Error::InvalidValue` and the existing curve is untouched. |
| `double eval(const double x) const` | `BSpline2d::eval(&self, f64) -> Result<f64>` | Returns `Ok(0.0)` when not `ok`, as the C++ returns 0. Never fails on range — see "Evaluation outside the domain". |
| `double derivative(const double x) const` | `BSpline2d::derivative(&self, f64) -> Result<f64>` | The eol-bspline name is `slope`. The fitted mean is deliberately not added, as in the source. |
| `bool ok() const` | `BSpline2d::ok(&self) -> bool` | Always true for a freshly constructed value, because construction failure is an `Err`. Retained because a failed `solve` can still make it false, and because the C++ contract "not ok, therefore evaluate to zero" is preserved. |
| `static void debug(bool enable)` | **not ported** | It sets a function-local `static bool` inside `BSplineBase<T>::Debug` that gates `std::cerr` tracing. The port emits no logging, and a process-global mutable flag is exactly the singleton the porting rules replace with caller-owned state. The class test marks the section `NOT_TESTABLE`. |
| `eol_bspline::BSpline<double>* spline_` (private) | private `Domain`, `Band` and coefficient vectors | Not public API. |

Native accessors, added because they expose internal state the PIMPL hides and
because the probe can check them:

| Rust | C++ counterpart |
|---|---|
| `BSpline2d::node_count()` | `BSplineBase::nNodes()` |
| `BSpline2d::node_spacing()` | `BSplineBase::DX` |
| `BSpline2d::domain()` | `(BSplineBase::Xmin(), BSplineBase::Xmax())`, the latter reconstructed as `Xmin + M*DX` exactly as the C++ does |
| `BSpline2d::data_range()` | the protected `xmin` / `xmax` members, which differ from `domain()` in the last bits |
| `BSpline2d::alpha()` | `BSplineBase::Alpha()` |
| `BSpline2d::wavelength()` | `BSplineBase::waveLength`, after the setup's rewrite |
| `BSpline2d::fitted_mean()` | `BSpline<T>::mean` |
| `BSpline2d::coefficient(n)` | `BSpline<T>::coefficient(int)`, including its "zero outside 0..=M or when not ok" contract |
| `BSpline2d::MAX_POINTS`, `MAX_NODES` | none; bounded-work ceilings |
| `impl SplineFunction for BSpline2d` | the duck typing `Math::spline_bisection` relies on |

## Preserved source conventions

These are the details that a "cleaner" reimplementation would have lost. Each is
reproduced deliberately.

1. **`PI` is `3.1415927`, not pi.** `BSplineBase.cpp:109` defines
   `const double BSplineBase<T>::PI = 3.1415927`. It enters
   `Alpha(wl) = ((wl / (2*PI*DX))^2)^2`, the weight of the whole derivative
   constraint. Using `std::f64::consts::PI` changes `alpha` in the eighth
   significant digit and every coefficient with it. The Rust constant is
   `EOL_PI` and carries an `#[allow(clippy::approx_constant)]` precisely because
   it is *not* pi.
2. **A zero cutoff wavelength does not disable the derivative constraint.** Both
   the OpenMS and the eol-bspline documentation say it does. `Setup` instead
   rewrites `waveLength = 1.0` — one unit of `x` — before `Alpha` is evaluated
   (`BSplineBase.cpp:560` for an explicit node count, `:566` for the automatic
   one). For the class-test fixture that leaves `alpha = 954.99`, not zero. The
   port reproduces the rewrite and contradicts the documentation at the item.
3. **The `float` accumulators.** `calculateQ` declares `float b1, b2, q` and
   `addP` declares `float pm, pn, sum`, so the boundary correction and every
   entry of the least-squares normal matrix are rounded to single precision
   before being accumulated into the `double` band. That is worth about seven
   significant digits and is visible in the coefficients. The Rust code rounds
   through `f32` at exactly the same points, including the `q += <double>`
   re-rounding after each addition.
4. **The `Setup` node search, including its starting point.** `ni` is seeded at
   9 and pre-incremented, so the first trial grid has ten intervals. The
   consequence, which no comment mentions: with an automatic node count and a
   positive cutoff wavelength, a dataset of ten or fewer points fails
   immediately, because `Ratiod` returns `NX/(ni+1) < 1`. The probe records
   `bspline_small_wl3 ok = 0` for eight points and a wavelength of 3. Also
   preserved: the second loop's use of `ratiof` from the *last* trial when
   deciding whether to continue, and its `--ni` before breaking.
5. **`Setup` only looks at the extremes of `x`.** The abscissae need not be
   sorted, and repeats are fine; only `min` and `max` set the grid. The class
   test's `ok()` section relies on this, fitting the same ramp forwards and
   backwards. The `else if` in the min/max scan is transcribed as written.
6. **Evaluation is a four-term sum over a clamped window.** `evaluate` computes
   `n = (int)((x - xmin)/DX)` and sums `A[i] * Basis(i, x)` for
   `i` in `max(0, n-1) ..= min(M, n+2)`, then adds the fitted mean. `slope` does
   the same with `DBasis` and does *not* add the mean.
7. **The fit is on mean-subtracted ordinates.** `BSpline<T>::solve` subtracts the
   mean of `y` before accumulating the right-hand side and `evaluate` adds it
   back, in that order and with that accumulation.
8. **The banded storage discards out-of-band writes.** `BandedMatrix::element`
   returns a shared `out_of_bounds` scratch reference for any coordinate outside
   the band or outside `[0, N)`, so a write there lands nowhere and a read
   returns whatever was written last. `Band` ignores such writes and reads them
   back as zero. The two are indistinguishable because no value the C++ leaves
   in that scratch element ever reaches an in-band coordinate: every
   `Q[i][j] += q` that lands outside also stores its result outside, and both
   `LU_factor_banded` and `LU_solve_banded` touch only `|i - j| <= 3` with both
   indices in `[1, N]`. For `N >= 3` the `Band::slot` test is literally the
   C++ `check_bounds` for a `setup(N, -3, 3)` matrix; the two-node grid of the
   next item is the one case where `N` is not `M + 1`.
9. **Two nodes leave the matrix at its default 1x1 shape, and the fit proceeds
   anyway.** `calculateQ` opens with `Q.setup(M+1, 3)`, which returns `false`
   for `M + 1 < 3` because the dimension is below the bandwidth, and the return
   value is ignored. `Q` therefore keeps the 1x1 shape its default constructor
   gave it. `calculateQ` still writes `qDelta(0, 0)` and both boundary
   corrections into that one element and `addP` still adds the normal-matrix
   entry to it, so it is not zero, `LU_factor_banded` succeeds on the 1x1
   system, `LU_solve_banded` divides the first right-hand-side entry by it, and
   `BSpline<T>::solve` sets `OK`. The second coefficient is never solved for:
   the right-hand side and the solution are the same `std::vector`, sized
   `M + 1 = 2`, while the matrix reports one row, so `A[1]` is left holding the
   raw `sum_j (y_j - mean) * Basis(1, x_j)` the solve accumulated.

   **The port reproduces this**, because it is what the source does and it is
   deterministic. For `x = 0..7`, `y = x^2`, `BC_ZERO_SECOND` and two nodes the
   probe records `ok = 1`, `alpha = 2.6723192544913614e-07`,
   `coeff0 = -8.0110109563847747`, `coeff1 = 15.000000000000011`,
   `eval(0) = 22.741741782711426`, `eval(3.5) = 23.711179967485478` and
   `derivative(3.5) = 2.0636797453269411`; `tests/spline_math.rs`
   (`bspline_two_node_grids_reproduce_the_degenerate_source_fit`) reproduces all
   of those bit-for-bit under three boundary conditions. It is not a fit — at
   3.5 the data is 12.25 and the second coefficient is in the units of the
   right-hand side rather than of the curve — and the type documentation says
   so, but a port does not get to decide that on the caller's behalf.

   Two caveats belong with those numbers. First, the C++ only reaches this state
   with `NDEBUG` set: on the two-node grid `Beta` is called outside the range its
   own `assert(0 <= m && m <= 3)` permits, and a build with assertions enabled
   aborts there instead of returning a spline. The oracle is a release build, as
   OpenMS ships. Second, the values do not depend on what those out-of-range
   reads returned, because each one is consumed by a write the banded storage
   discards — which is also why the port can skip the reads and still agree.

   Three nodes is the smallest grid whose matrix the C++ actually shapes
   (`setup(3, -3, 3)` is accepted, its outermost bands having length zero), and
   it is a real fit: `eval(3.5) = 11.676270573337504` for the same data.
10. **The LU code is transcribed one-based.** `lu_factor_banded` and
    `lu_solve_banded` keep the source's `A(i, j)` indexing, its band-limited
    inner loops and its accumulation order, including the unchecked
    `b[M-1] /= A(M, M)` that precedes the checked back-substitution loop.

## Native differences

| Difference | Source behaviour | Port behaviour | Why |
|---|---|---|---|
| Failed setup or factoring | Object exists, `ok()` is false, `eval` returns 0 for every `x`, and reading `nNodes()` or `Alpha()` returns uninitialised members. The probe shows exactly that for `bspline_sinus_bc2_wl100` and `bspline_small_wl3`. | `Err(Error::InvalidValue)` from the constructor, with a message naming the condition. | An object whose every accessor is meaningless is what `Result` is for; and the C++ reads uninitialised `M` and `DX`, which is undefined behaviour a caller cannot detect. |
| `solve` failure | `OK = false`, coefficients left as partial garbage, `eval` returns 0. | `Err`, coefficients zeroed, `ok()` false, `eval` returns `Ok(0.0)`. The C++ "evaluate to zero" contract is kept. | Atomicity: the new curve is built in a temporary and only committed on success, so a rejected `solve` — for instance one with the wrong length — leaves the previous curve exactly as it was. |
| `solve` length mismatch | Inert precondition, then reads `NX` elements from a shorter array. | `Err` before anything is touched. | Out-of-bounds read. |
| Non-finite `x` in `eval` / `derivative` | `(int)((x - xmin)/DX)` on a `NaN` is undefined behaviour; in practice it yields node 0 and the function returns the fitted mean. | `Err(Error::InvalidValue)`. | Undefined behaviour, and a mean returned for a `NaN` query is worse than an error. |
| Non-finite coefficients | `LU_solve_banded` can divide by a zero pivot at `b[M-1] /= A(M,M)`, which it does not check; `OK` is then set true and every evaluation is `NaN`. | The solution is checked for finiteness; a non-finite coefficient is an `Err` and leaves the spline not `ok`. | No silent `NaN` propagation. |
| Out-of-range `Beta` | `calculateQ`'s lower-right loop calls `Beta(M-4)`, negative when `M < 4`, and its upper-left loop calls `Beta(j)` for `j` up to `4`, out of range whenever `M <= 3`. With `NDEBUG` both read outside `BoundaryConditions`; with assertions enabled they abort. | Those `(i, j)` pairs are skipped. | **Numerically identical**, and this is the load-bearing part of the argument: in every such case the pair also lies outside the banded matrix (`min(i, j) < 0`, or `max(i, j) > M`, or — on the two-node grid — outside the 1x1 shape the failed `setup` left behind), so both the `Q[i][j] += q` and the `Q[j][i] = ...` that follow write into the discarded scratch element. The correction is computed and thrown away. Skipping it removes the out-of-bounds read and changes nothing. `BSplineSmoothingSpline` reaches this on every dataset, because its candidate list always tries four nodes. |
| Point and node counts | Unbounded. | `MAX_POINTS = 250 000`, `MAX_NODES = 500 001`, checked before allocation. The automatic node count for a wavelength-free fit is `2n + 1`, so the two ceilings are consistent. | Bounded work: the banded matrix is seven `f64` per node, so the cap is 28 MiB. |
| Degenerate abscissae | All `x` equal gives `DX = 0` and `NaN` everywhere. | `Err`. | Division by zero. |
| `debug(bool)` | Process-global mutable flag gating `std::cerr`. | Not ported. | See the API table. |

No OpenMP: neither `BSpline2d.cpp` nor any eol-bspline file carries a
`#pragma omp`, so there is no parallelism gap to record. The port is serial and
so is the source.

## Evaluation outside the fitted domain

This is the behaviour most likely to surprise a caller, so it is stated
precisely. `eval` never fails on range. Each basis function is zero more than
two node intervals from its node, and the evaluation window is clamped to
`[0, M]`; beyond `xmin - 2*DX` and `xmax + 2*DX` the window is empty and the sum
is just the fitted mean. Between the domain edge and that point the curve decays
smoothly towards the mean. `derivative` likewise returns exactly zero out there.

The probe pins it: for `y = x^2` on `x = 0..7`, `eval(-2) == eval(9) == 17.5`,
which is `140/8`, the mean of the ordinates — a value the Rust test derives
arithmetically rather than transcribing.

`CubicSpline2d` does the opposite and rejects any query outside its knot range.
A caller switching between the two must handle that difference explicitly.

## Checked boundaries and evidence

Evidence tier 2 (executed probe), bit-for-bit. The probe compiles
`BSpline2d.cpp` and the six vendored eol-bspline files unmodified against shims
for the build-configuration headers, and additionally instantiates
`eol_bspline::BSpline<double>` directly so the hidden domain and coefficients
are observable. Flags, compiler, platform and sha256 values are in
`tests/data/spline_math_provenance.json`.

**Floating-point contraction.** The probe is built with `-ffp-contract=off`.
Two strict builds at `-O0` and `-O2` produce byte-identical output, so the
recorded oracle is optimisation-independent, and the port matches it exactly.

Built with clang's default contraction on arm64, which fuses multiply-adds into
FMA, the same sources produce a different number in 1 234 of the 1 828 recorded
values. An earlier version of this document called that difference "one or two
units in the last place"; it is not, and the measured figures are these. The
median relative difference is `8.9e-15`, 300 values differ by more than `1e-13`
relative and 35 by more than `1e-12`; the largest gaps on values of ordinary
magnitude are `9.7e-9` relative on `bspline_ramp_asc eval@0` and `1.7e-10` on
`bspline_solve eval@0`. Fourteen rows differ by order unity, and every one of
them is a quantity that is numerically zero — the residual sum of squares of an
exactly fitted line or parabola, or its value at a point where the curve crosses
zero; for instance `smooth_linear5_auto rss` is `0` strict and `5.1e-31`
contracted. The banded normal equations are the amplifier: the ill-conditioned
solve turns a fused rounding into a coefficient difference in the fourteenth
digit, and `BSplineSmoothingSpline` inherits it.

That is a property of the C++, recorded here rather than hidden behind a
tolerance. The consequence for a consumer is that the port agrees bit-for-bit
with a strict C++ build and to about thirteen significant digits with a
contracted one — and that a residual sum of squares near zero should not be
compared exactly across the two.

### Class-test sections

`BSpline2d_test.cpp` has 7 `START_SECTION`s and 8 assertion macros (one of them
`NOT_TESTABLE`). All 7 are mapped.

| Section | Rust test | One value it reproduces |
|---|---|---|
| `BSpline2d(const std::vector<double>&, const std::vector<double>&, double, BoundaryCondition, Size)` | `bspline_reproduces_the_upstream_sinus_fixture` | The default fit of the 202-point fixture has `node_count() == 405` and `coefficient(0) == 2.6003560634855045`; all 405 coefficients are compared. The section's fourth construction, `BSpline2d(x, y, 100, BC_ZERO_SECOND)`, cannot succeed — 100 exceeds the 11.57 span of the abscissae — and the same test asserts the probe's `ok = 0` alongside the port's `Err`. |
| `virtual ~BSpline2d()` | not a Rust construct | The C++ destructor only deletes the PIMPL pointer. The Rust type owns `Vec<f64>` and has no `Drop` impl; `bspline_solve_refits_the_same_domain` and the other tests construct and drop many splines. Mapped as a native equivalent with no behaviour to assert. |
| `bool solve(const std::vector<double>&)` | `bspline_solve_refits_the_same_domain` | After refitting the fixture to `10*sin(x)`, `eval(-8.0) == 1.7716867914263408`. |
| `double eval(const double x) const` | `bspline_smoothing_beats_the_noise_by_the_margin_the_class_test_requires` | Mean squared error against `10*sin(x)` falls from `1.3652888205904805` to `0.18800316671802897` for the default fit and `0.13867348442977684` at wavelength 2, and the section's own assertion — smoothed error below half the noisy error — is reproduced. |
| `double derivative(const double x) const` | the same test | Mean absolute derivative error `1.0649541454021494`, and the section's `< 2.0` assertion. |
| `bool ok() const` | `bspline_is_insensitive_to_the_order_of_the_abscissae` | The ascending ramp gives `eval(5.5) == 5.4999999995195745` and the descending one `5.4999999999999991`; the section's `TEST_REAL_SIMILAR(y1, y2)` holds, and all 100 coefficients of each fit are compared. |
| `void debug(bool enable)` | not ported | The section is `NOT_TESTABLE`; see the API table for why the static flag is not carried across. |

Unaccounted sections: none.

### Beyond the class test

The probe covers seventeen configurations, all compared coefficient by
coefficient: the fixture at five wavelength and boundary-condition combinations,
an eight-point quadratic at automatic, four- and six-node grids under two
boundary conditions, a 21-point quadratic at two cutoff wavelengths, the ramp
forwards and backwards at 100 nodes, the same eight-point quadratic on the
degenerate two-node grid under all three boundary conditions and on the
three-node grid, and two setups that fail.

### Independently derived checks

* `eval` far outside the domain equals the arithmetic mean of the ordinates, and
  `derivative` there is exactly zero.
* Under `BC_ZERO_FIRST` the slope at both ends of the node domain is zero — `0.0`
  exactly at the lower end and `2.3e-13` at the upper — which the test asserts as
  a property of the boundary condition, not as a transcribed number.
* Negating every ordinate negates the fitted mean and every value, because the
  normal equations are linear in `y`.
* A rejected `solve` leaves the previous curve bit-identical.
* On the two-node grid the unsolved second coefficient equals the right-hand
  side the solve accumulated, recomputed in the test from the basis and the
  mean-subtracted ordinates rather than transcribed, and negating the ordinates
  negates it.

### Boundaries the port checks

* `x` and `y` the same length, non-empty, all finite, at most `MAX_POINTS`.
* `wavelength` finite and non-negative; `num_nodes` at most `MAX_NODES`.
* At least two distinct abscissae, so the node spacing is positive and finite.
* The derived node count at most `MAX_NODES`, with the search loop bounded by
  the same ceiling. There is no lower bound on the node count: two nodes is the
  degenerate grid of item 9 and is accepted, exactly as the source accepts it.
* The penalty weight finite.
* The banded factorisation succeeding, and the solution finite.
* Query positions finite; results finite.
