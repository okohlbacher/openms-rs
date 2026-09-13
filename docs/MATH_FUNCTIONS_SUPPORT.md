# General math functions

`concept::math_functions` ports the `OpenMS::Math` namespace of
`src/openms/include/OpenMS/MATH/MathFunctions.h` at revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`. The header is header-only: its
paired `src/openms/source/MATH/MathFunctions.cpp` opens and closes
`namespace OpenMS` and defines nothing, so every inline body in the header *is*
the implementation and was read as such.

`MathFunctions.h` has eleven direct TOPP consumers - `DecoyDatabase`,
`FeatureFinderMetaboIdent`, `FileInfo`, `HighResPrecursorMassCorrector`,
`IDExtractor`, `MassTraceExtractor`, `MzTabExporter`, `OpenNuXL`,
`OpenSwathAssayGenerator`, `QCMerger` and `TextExporter` - the highest fan-out in
the MATH domain. The two families those consumers reach for are the ppm/Dalton
tolerance conversions and the decimal rounding helpers, and both are where the
port is most conservative.

Rust files covering the header: `src/concept/math_functions.rs` (one module).
Tests: `tests/math_functions.rs` plus three private unit tests in the module.

## Why `concept` and not `math`

The header's own Doxygen puts the `Math` namespace in `@ingroup Concept`, and
`src/math/` does not exist in this crate. Creating it would mean adding a
top-level `pub mod math;` to `src/lib.rs`, which this work package may not edit;
`src/concept/` is the closest existing home and the one the source's own grouping
names. The module reaches for nothing but `crate::{Error, Result}`, so
`tools/check_module_cycles.py` reports the same 57 edges and 14 mutually
dependent pairs as before.

## API mapping

Every public member of the header appears below.

| C++ member | Rust counterpart | Notes |
|---|---|---|
| `template<T> bool extendRange(T& min, T& max, const T& value)` | `extend_range(&mut f64, &mut f64, f64) -> bool` | In/out references kept; the return value alone cannot say which bound moved. `f64` only. |
| `template<T> bool contains(T value, T min, T max)` | `contains(f64, f64, f64) -> bool` | `min <= value && value <= max`, NaN answers `false`. |
| `pair<double,double> zoomIn(double left, double right, float factor, float align)` | `zoom_in(f64, f64, f32, f32) -> Result<(f64, f64)>` | `factor`/`align` stay `f32` because `(1.0f - factor)` is a single-precision subtraction in the source. Debug-only preconditions become checked errors. |
| `using BinContainer = std::vector<RangeBase>` | `Vec<Bin>` | See `Bin` below. |
| `BinContainer createBins(double min, double max, uint32_t number_of_bins, double extend_margin = 0)` | `create_bins(f64, f64, u32, f64) -> Result<Vec<Bin>>` | Default argument becomes an explicit `0.0`. `MAX_BINS` ceiling added. |
| `double ceilDecimal(double x, int decPow)` | `ceil_decimal(f64, i32) -> Result<f64>` | Same `pow(10.0, decPow)` scaling. |
| `double roundDecimal(double x, int decPow)` | `round_decimal(f64, i32) -> Result<f64>` | Both branches transcribed, including the negative-zero result at `x == 0`. |
| `double intervalTransformation(double x, double left1, double right1, double left2, double right2)` | `interval_transformation(f64, f64, f64, f64, f64) -> Result<f64>` | Multiplication-before-division order preserved; zero source span checked. |
| `double linear2log(double x)` | `linear_to_log10(f64) -> Result<f64>` | Renamed for Rust; `x + 1 <= 0` checked. |
| `double log2linear(double x)` | `log10_to_linear(f64) -> Result<f64>` | Renamed for Rust. |
| `bool isOdd(UInt x)` | `is_odd(u32) -> bool` | `UInt` is `u32`. |
| `template<T> T round(T x)` | `round(f64) -> f64` | `std::round` and `f64::round` share the away-from-zero tie rule. The class test's `float` instantiation is `f32::round`, asserted in the test but not wrapped. |
| `template<T> T roundTo(const T value, int digits)` | `round_to(f64, i32) -> Result<f64>` | The source's factor **loop** is transcribed, not replaced by `powi`. `MAX_DECIMAL_DIGITS` ceiling added. |
| `template<T> double percentOf(T value, T total, int digits)` | `percent_of(f64, f64, i32) -> Result<f64>` | `Exception::InvalidValue` becomes `Error::InvalidValue`; the `total <= 0` shortcut is kept. |
| `bool approximatelyEqual(double a, double b, double tol)` | `approximately_equal(f64, f64, f64) -> bool` | Absolute, not relative. |
| `template<T> T gcd(T a, T b)` | `gcd(i64, i64) -> Result<i64>` | `i64` only. `checked_rem` replaces the source's undefined `INT_MIN % -1`. |
| `template<T> T gcd(T a, T b, T& u1, T& u2)` | `extended_gcd(i64, i64) -> Result<ExtendedGcd>` | Out-parameters become the `ExtendedGcd { gcd, u1, u2 }` return. All arithmetic checked. |
| `template<T> T getPPM(T mz_obs, T mz_ref)` | `ppm(f64, f64) -> Result<f64>` | Signed; divides by the **reference**. Zero reference checked. |
| `template<T> T getPPMAbs(T mz_obs, T mz_ref)` | `ppm_abs(f64, f64) -> Result<f64>` | `|ppm|`. |
| `template<T> T ppmToMass(T ppm, T mz_ref)` | `ppm_to_mass(f64, f64) -> Result<f64>` | Signed; `(ppm / 1e6) * mz_ref` in that order. |
| `template<T> T ppmToMassAbs(T ppm, T mz_ref)` | `ppm_to_mass_abs(f64, f64) -> Result<f64>` | Its C++ comment block opens with `/*`, so this member is absent from the generated Doxygen; its content is carried into the rustdoc anyway. |
| `pair<double,double> getTolWindow(double val, double tol, bool ppm)` | `tolerance_window(f64, f64, bool) -> Result<(f64, f64)>` | The asymmetric ppm window is reproduced exactly; `tol == 1e6` is checked instead of dividing by zero. |
| `template<T1> T1::value_type quantile(const T1& x, double q)` | `quantile(&[f64], f64) -> Result<f64>` | `f64` samples only. `Exception::InvalidParameter` on empty becomes `Error::InvalidValue`; sortedness is now checked (`Error::UnsortedData`). |
| `class RandomShuffler` | **not ported here**: `chemistry::decoy_random::DecoyRandom` already implements exactly these semantics (the `boost::mt19937_64` word stream, taken from `rand_mt::Mt64`; Boost `uniform_int` bucket mapping; descending Fisher-Yates) but is `pub(crate)` inside another group's module and specialised to `&mut [u8]`. See "Deferred: RandomShuffler" below. |
| `RandomShuffler::RandomShuffler(int seed)` | not ported here; `DecoyRandom::seeded(u64)` | |
| `RandomShuffler::RandomShuffler(const boost::mt19937_64&)` | not ported: constructing from a foreign engine value has no Rust meaning without that engine type. | |
| `RandomShuffler::RandomShuffler()` / `~RandomShuffler()` | not ported: a default engine state and a trivial destructor have no counterpart. | |
| `RandomShuffler::rng_` (public member) | not ported: the engine state stays private. Exposing a mutable engine field is exactly the "global mutable singleton" shape the port avoids. | |
| `template<RandomAccessIterator> void RandomShuffler::portable_random_shuffle(first, last)` | not ported here; `DecoyRandom::shuffle(&mut [u8], &mut usize)` | The generic element type is the part that is missing. |
| `void RandomShuffler::seed(uint64_t val)` | not ported here; `DecoyRandom::reseed(u64)` | |
| `double log_binomial_coef(unsigned n, unsigned k)` | `log_binomial_coef(u32, u32) -> Result<f64>` | Both shortcuts and the `k > n/2` integer-division swap kept. `boost::math::lgamma` becomes `libm::lgamma`. |
| `double log_sum_exp(double x, double y)` | `log_sum_exp(f64, f64) -> f64` | Total, as in the source; `std::max` semantics reproduced rather than `f64::max`. |
| `double binomial_cdf_complement(unsigned N, unsigned n, double p)` | `binomial_cdf_complement(u32, u32, f64) -> Result<f64>` | All four shortcuts kept; the tail itself uses a different algorithm, stated below. |

Native additions, with no source counterpart:

| Rust item | Purpose |
|---|---|
| `Bin { min, max }`, `Bin::is_empty`, `Bin::contains` | The two coordinates of `RangeBase` that `createBins` actually uses, plus the two predicates its body depends on. `kernel::ranges::RangeBase` is not reachable from `concept` without a new module-pair edge. |
| `ExtendedGcd { gcd, u1, u2 }` | Return type replacing two out-parameters. |
| `MAX_BINS`, `MAX_BINOMIAL_TRIALS`, `MAX_QUANTILE_ITEMS`, `MAX_DECIMAL_DIGITS` | Explicit ceilings on the four operations whose cost scales with a caller-supplied number. |

## Preserved source conventions

These are the places where a "cleaner" rewrite would change results, so the
source spelling is kept verbatim.

* **ppm is anchored on the reference mass.** `getPPM` divides by `mz_ref`, never
  by `mz_obs`, so the function is not antisymmetric: `ppm(1000, 1001)` is
  `-999.000999…`, a full part per million away from the `-1000` a symmetric
  definition would give. Swapping the two arguments at a call site is therefore
  a silent, mass-dependent error rather than a sign flip, and
  `ppm_divides_by_the_reference_and_is_not_antisymmetric` asserts exactly that.
* **The ppm tolerance window is asymmetric.** `getTolWindow` computes
  `left = val - val*tol*1e-6` but `right = val / (1 - tol*1e-6)`. The right edge
  is the largest `x` that still has `val` inside `x`'s own ppm window, which is
  what makes the compatibility relation between two masses symmetric even though
  the interval is not. `tolerance_window_is_asymmetric_in_ppm_mode_and_symmetric_in_dalton_mode`
  asserts `right - val > val - left` and confirms both edges are exactly at the
  requested ppm from the other endpoint.
* **`ppmToMass` divides before it multiplies.** `(ppm / 1e6) * mz_ref` is not the
  same double as `ppm * mz_ref / 1e6`, so the order is preserved.
* **`roundTo` builds its factor in a loop.** The source multiplies or divides a
  running `1.0` by ten `|digits|` times instead of calling `pow`. That is a
  different double from `10f64.powi(digits)` for most magnitudes, and
  `round_to_uses_the_source_scaling_loop` pins the port to the loop value bit for
  bit while the `powi` form only agrees within the class test's tolerance.
* **`ceilDecimal` and `roundDecimal` use `pow(10.0, decPow)` with a promoted
  `int`.** That selects `std::pow(double, double)`, so `f64::powf` is the exact
  mapping; `powi` evaluates a different multiplication sequence.
* **`roundDecimal(0.0, p)` returns negative zero.** Zero is not greater than
  zero, so it takes the negating branch. Reproduced rather than normalised,
  because `copysign` and division make it observable.
* **`log_binomial_coef`'s symmetry swap.** `k > n / 2` uses integer division and
  replaces `k` with `n - k` *before* the log-gammas, which makes `C(10,3)` and
  `C(10,7)` bit-identical rather than merely close. The test asserts exact
  equality.
* **`log_sum_exp` uses `std::max`, not `f64::max`.** `std::max(a, b)` is
  `(a < b) ? b : a` and returns its first argument when the comparison is false,
  so a NaN `x` propagates; `f64::max` would discard it. The private unit test
  `log_sum_exp_uses_the_source_maximum_not_the_rust_one` pins this.
* **`quantile`'s index convention.** The source index is `max(0, n*q - 1)`, not
  the usual `(n - 1) * q`. For `[1,2,3,4,5]` and `q = 0.5` it interpolates
  between `x[1]` and `x[2]` and returns `2.5`, not `3`. Reproduced, because
  callers are calibrated against it, and called out here and in the rustdoc
  because it is easy to mistake for a textbook quantile.
* **`binomial_cdf_complement`'s four shortcuts.** `n == 0` returns `1.0` before
  anything else, `p == 0` returns `0.0` for a positive `n`, `p == 1` returns
  `1.0`. These are exact returns, not limits of the general formula.
* **`createBins` repairs the outer borders after extending.** `setMin`/`setMax`
  also move the opposite bound when the reset would invert the bin; that repair
  is transcribed.

## Native differences

Every divergence, and why.

1. **Checked instead of unchecked arithmetic.** The source divides by
   `mz_ref`, by `right1 - left1`, by `1 - tol*1e-6` and by the rounding factor
   without testing any of them, and calls `log10(x + 1)` for any `x`. Each is
   checked here and returns `Error::InvalidValue` rather than an infinity or a
   NaN, because these results flow straight into mass matching where a NaN
   compares false against every tolerance and silently drops a match.
2. **Non-finite inputs are refused** wherever a result would otherwise be a NaN.
   The one exception is `log_sum_exp`, which is total in the source and stays
   total here so that the `-inf` identity remains usable for log-domain
   accumulation.
3. **Preconditions become errors.** `zoomIn` and `createBins` state their
   preconditions with `OPENMS_PRECONDITION`, which is compiled out of release
   builds. This port checks unconditionally: `Error::InvalidRange` for
   `min >= max`, `Error::InvalidValue` for a zero bin count, a negative factor or
   an out-of-range alignment.
4. **`quantile` checks sortedness** and returns `Error::UnsortedData`. The source
   states the precondition in its `@brief` only. The cost is one pass and the
   alternative is a number that is not a quantile of anything.
5. **Integer overflow is checked** in both `gcd` overloads. `gcd(i64::MIN, -1)`
   is undefined behaviour in C++ (and traps on x86); here it is an error.
6. **`binomial_cdf_complement` uses a different algorithm**, and this is the one
   place where the port does not reproduce the source's evaluation. The source
   calls `boost::math::cdf(complement(binomial_distribution<double>(N, p), n-1))`,
   which is the regularized incomplete beta function evaluated by a continued
   fraction. Boost is not a dependency of this crate and no dependency may be
   added here, so the tail is summed as its own definition instead:
   `sum over k in n..=N of C(N,k) p^k (1-p)^(N-k)`, evaluated in the log domain
   with the largest term factored out (`exp(max) * sum exp(t_k - max)`) and
   accumulated in ascending `k`. `ln(1-p)` uses `ln_1p(-p)` for small `p`.
   *Agreement:* against an exact rational evaluation of the same sum (Python
   `fractions.Fraction` with `math.comb`, an independent oracle, not a C++ run)
   the largest relative deviation over all nine numeric cases in the source's
   class test is `2.5e-14`, at `N = 100`; the dyadic case `B(10, 0.5)` agrees
   with `638/1024` to `3.4e-15`. Both are four to nine orders of magnitude
   inside the class test's own `1e-5` relative tolerance. The cost is linear in
   `N` rather than near-constant, which is why `MAX_BINOMIAL_TRIALS` exists.
   The result is clamped into `[0, 1]`, because the factored sum can overshoot
   one by a few units in the last place while a probability cannot.
7. **`libm::lgamma` replaces `boost::math::lgamma`** in `log_binomial_coef`.
   These are different implementations of the same function; over this domain
   they agree to a few units in the last place. `ln(C(10,5)) = ln(252)` and
   `ln(C(20,10)) = ln(184756)` are asserted to `1e-14` relative against the
   closed form, which is an independent check rather than a transcribed literal.
8. **Template parameters are fixed.** `extendRange`, `contains`, `round`,
   `roundTo`, `percentOf`, `getPPM*`, `ppmToMass*` and `quantile` are templates;
   the port instantiates `f64` (and `i64` for `gcd`), since those are the only
   instantiations the class test and the TOPP consumers use. `round`'s `float`
   instantiation is `f32::round` in Rust and needs no wrapper; the test asserts
   both source float literals against it.
9. **`createBins` returns `Bin`, not `RangeBase`.** See the API table.
10. **Bounded work.** Four ceilings are added where the source has none:
    `MAX_BINS` (before the vector is allocated), `MAX_BINOMIAL_TRIALS` (before the
    term vector is allocated), `MAX_QUANTILE_ITEMS` (before the ordering scan) and
    `MAX_DECIMAL_DIGITS` (before `roundTo`'s factor loop, which the source will
    happily run two billion times for `roundTo(x, INT_MIN)`).
11. **No OpenMP.** The header has no `#pragma omp`, so there is no parallelism
    gap to record for this file.

## Deferred: RandomShuffler

`Math::RandomShuffler` is the only public member of the header without a
counterpart in this module, and the reason is placement rather than difficulty.
Its exact semantics are already implemented, with an independent oracle, in
`src/chemistry/decoy_random.rs` as `DecoyRandom`. The `boost::mt19937_64` word
stream comes from the `rand_mt::Mt64` dependency; Boost `uniform_int`'s
bucket-division-with-rejection mapping and the *descending* Fisher-Yates loop
are ported beside it. `Mt64` omits Boost's seed-time state normalization, which
changes no output word; `DECOY_REFERENCE_REVIEW.md` gives the reason and the
executed comparison against the hand-written engine it replaced. That type is
`pub(crate)`, its `shuffle` is `pub(super)`, and it shuffles `&mut [u8]`.

Publishing a generic shuffler from `concept::math_functions` would mean one of
two things, and neither belongs in this work package:

* duplicating the Boost range mapping and the Fisher-Yates loop, which would put
  two copies of bit-exactness-critical code in the crate (the engine is a shared
  dependency and needs no copy), or
* widening the visibility of another group's private module and generalising its
  element type, which changes `DECOY_SUPPORT`'s documented surface and its
  ledger scope.

So it is reported as a deferral for the integrator, with the recommendation to
promote `DecoyRandom` into a shared module in a package that owns both files.
Until then, `chemistry::decoy_generator` is the only consumer, and it is the only
consumer upstream too apart from `DecoyDatabase`.

## Checked boundaries and evidence

`MathFunctions_test.cpp` has **16** `START_SECTION` blocks and all 16 are mapped;
the table below cites the Rust test and one concrete value it reproduces.

| Section | Rust test | Cited value |
|---|---|---|
| `log_binomial_coef` | `log_binomial_coef_matches_source_values_and_is_exactly_symmetric` | `log_binomial_coef(10, 5) == 5.5294` |
| `log_sum_exp` | `log_sum_exp_is_stable_and_treats_negative_infinity_as_an_identity` | `log_sum_exp(1.0, 2.0) == 2.31326169` |
| `ceilDecimal` | `ceil_decimal_shifts_ceils_and_shifts_back` | `ceil_decimal(12345.67, 1) == 12350.0` |
| `roundDecimal` | `round_decimal_rounds_to_nearest_and_keeps_the_negative_zero_branch` | `round_decimal(12345.67, 2) == 12300.0` |
| `intervalTransformation` | `interval_transformation_maps_between_two_spans` | `interval_transformation(0.5, 0.25, 1.0, 0.0, 600.0) == 200.0` |
| `linear2log` | `linear_to_log10_adds_one_before_taking_the_logarithm` | `linear_to_log10(99.0) == 2.0` |
| `log2linear` | `log10_to_linear_inverts_linear_to_log10` | `log10_to_linear(3.0) == 999.0` |
| `isOdd` | `is_odd_tests_the_low_bit` | `is_odd(3) == true` |
| `round` | `round_breaks_ties_away_from_zero` | `round(-675.77) == -676.0` |
| `roundTo` | `round_to_uses_the_source_scaling_loop` | `round_to(1234.9, -2) == 1200.0` |
| `percentOf` | `percent_of_rounds_and_refuses_negative_arguments` | `percent_of(1/3, 1.0, 4) == 33.3333` |
| `approximatelyEqual` | `approximately_equal_is_an_absolute_comparison` | `approximately_equal(1.1, 1.1002, 0.0001) == false` |
| `getPPM` | `ppm_divides_by_the_reference_and_is_not_antisymmetric` | `ppm(999.0, 1000.0) == -1000.0` |
| `getPPMAbs` | `ppm_abs_drops_the_sign_only` | `ppm_abs(999.0, 1000.0) == 1000.0` |
| `getTolWindow` | `tolerance_window_is_asymmetric_in_ppm_mode_and_symmetric_in_dalton_mode` | `tolerance_window(500, 5, true).1 == 500.0025000125` |
| `binomial_cdf_complement` | `binomial_cdf_complement_matches_source_values_and_an_exact_rational_oracle` | `binomial_cdf_complement(100, 60, 0.5) == 0.02844` |

Seven public members have **no** upstream section: `extendRange`, `contains`,
`zoomIn`, `createBins`, both `gcd` overloads, `ppmToMass`, `ppmToMassAbs`,
`quantile` and `RandomShuffler`. Their tests
(`extend_range_and_contains_follow_the_source_comparison_order`,
`zoom_in_round_trips_and_keeps_the_single_precision_subtraction`,
`create_bins_partitions_overlaps_and_never_extends_the_outer_borders`,
`gcd_and_extended_gcd_follow_knuth_and_check_their_arithmetic`,
`ppm_to_mass_carries_the_sign_and_inverts_the_ppm_computation`,
`quantile_uses_the_source_index_convention_and_checks_its_precondition`) are
independently derived from the header's documented behaviour, its inline bodies
and its own `@code` examples - which is tier 4 evidence, stronger than a
transcribed literal, but it means upstream has no assertion for any of them.

Evidence tier for the transcribed expected values is **3 (source review)**: no
C++ was built or executed, and none of the ten prebuilt reference outputs this
repository retains exercises `MathFunctions.h` directly. Where a value could be
derived independently it was, and the test says so inline:

* `binomial_cdf_complement` against exact rational arithmetic (tier 4).
* `log_binomial_coef` against `ln(252)` and `ln(184756)` (tier 4).
* `log_sum_exp(x, x) == x + ln 2` exactly (tier 4).
* The ppm round trip `ppm(ref + ppm_to_mass(t, ref), ref) == t` (tier 4).
* The `zoomIn` round trip from the header's own `@code` block (tier 4).
* `create_bins` adjacency, union and `2 * margin` overlap invariants (tier 4).
* The Bezout identity `a*u1 + b*u2 == gcd` for five pairs (tier 4).
* Endpoint invariance of `interval_transformation` (tier 4).

Resource behaviour: every function is `O(1)` except `create_bins` (`O(bins)`,
one allocation, preflighted against `MAX_BINS`), `binomial_cdf_complement`
(`O(trials - successes)` time and one allocation of the same size, preflighted
against `MAX_BINOMIAL_TRIALS`), `quantile` (`O(n)` scan, no allocation,
preflighted against `MAX_QUANTILE_ITEMS`) and `round_to`/`percent_of`
(`O(|digits|)`, no allocation, preflighted against `MAX_DECIMAL_DIGITS`). No
function mutates a caller's data except `extend_range`, which is the source's own
in/out contract and which on a NaN leaves both bounds untouched.

## Issues found in the C++

Recorded here for the integrator; this package may not edit
`OpenMS_CPP_ISSUES.md`.

1. **`createBins(min, max, 0, …)` is undefined behaviour in release builds.**
   `OPENMS_PRECONDITION(number_of_bins >= 1, …)` is inactive there, `res` is then
   an empty vector, and `res.front()` / `res.back()` are called on it
   unconditionally (`MathFunctions.h:129-137`). The division `(max-min)/0` is
   merely infinite; the `front()` is the actual defect.
2. **`Math::quantile` is off by one relative to the usual convention.** Its index
   is `max(0., n*q - 1)` (`MathFunctions.h:467`), so the median of a five-element
   sample is the midpoint of the second and third values rather than the third.
   Whether this is intended is not recorded anywhere in the header.
3. **`gcd(T a, T b)` and `gcd(T, T, T&, T&)` have undefined behaviour at the
   signed extremes** (`MathFunctions.h:316`, `:348`): `a % b` and `u3 / v3` with
   `INT_MIN` and `-1` overflow, and on x86 the hardware traps rather than
   wrapping.
4. **`binomial_cdf_complement` lets a NaN `p` through.** `p < 0.0 || p > 1.0` is
   false for a NaN (`MathFunctions.h:570`), so the value reaches Boost's
   distribution constructor, whose behaviour then depends on the configured
   policy rather than on the documented `std::invalid_argument`.
5. **`log_sum_exp(+inf, y)` returns NaN.** Only the negative infinities are
   guarded (`MathFunctions.h:548-549`); `max_val - max_val` is then `inf - inf`.
6. **Two members throw `std::invalid_argument`, not an OpenMS exception.**
   `log_binomial_coef` and `binomial_cdf_complement` (`MathFunctions.h:521`,
   `:572`) are the only members of the namespace that do, so a caller catching
   `Exception::BaseException` misses them.
7. **`ppmToMassAbs`'s documentation block opens with `/*` instead of `/**`**
   (`MathFunctions.h:410`), so the member is missing from generated Doxygen.
