# ZhangSimilarityScore support

Port of `src/openms/include/OpenMS/COMPARISON/ZhangSimilarityScore.h` and
`src/openms/source/COMPARISON/ZhangSimilarityScore.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `b55d2021fc4cea74dea52e27c97e2f8ca2a3620b856f374aeeaa0f594045fa12`,
`.cpp` sha256 `66aa29df1e5a4ce270a8b2e69aa2849b190fe6071c71f2c4a05b70183929c7d3`).

Rust: [`comparison::ZhangSimilarityScorer`](../src/comparison.rs), implementing
[`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class ZhangSimilarityScore : public PeakSpectrumCompareFunctor` | `pub struct ZhangSimilarityScorer` implementing `PeakSpectrumCompareFunctor` | composition instead of inheritance |
| `ZhangSimilarityScore()` | `ZhangSimilarityScorer::new() -> Result<Self>` | reproduces `ZhangSimilarityScore.cpp:20-32`: base handler, `setName("ZhangSimilarityScore")` at `:23`, four defaults, `defaultsToParam_()` |
| `ZhangSimilarityScore(const ZhangSimilarityScore& source)` | `Clone` | C++ copy constructor is `= default` |
| `~ZhangSimilarityScore() override` | drop glue | C++ destructor is `= default` |
| `ZhangSimilarityScore& operator=(const ZhangSimilarityScore& source)` | assignment of a clone | |
| `double operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const override` | `PeakSpectrumCompareFunctor::score` | `Result<f64>` |
| `double operator()(const PeakSpectrum& spec) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override at `:47-50` is `return operator()(spec, spec)` |
| `protected: double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const` | `ZhangSimilarityScorer::factor(mz_tolerance, mz_difference, is_gaussian) -> Result<f64>` | promoted from protected to a public associated function: it is a pure function of its three arguments, the port has no subclasses to protect it for, and exposing it lets the weighting be tested directly. Returns `Result` because a non-positive tolerance divides by zero upstream. The C++ default argument `is_gaussian = false` has no Rust equivalent; the only call site passes it explicitly |
| inherited `getParameters` / `setParameters` / `getName` | `handler()`, `handler_mut()`, `name()` | |
| - | `ZhangSimilarityScorer::max_pairs` | native resource ceiling, default [`DEFAULT_SCORED_PAIRS`](../src/comparison.rs) = 5000000 |
| `@htmlinclude OpenMS_ZhangSimilarityScore.parameters` | the parameter table in the rustdoc | |

Parameters: `tolerance` `0.2`, `is_relative_tolerance` `"false"`,
`use_linear_factor` `"false"`, `use_gaussian_factor` `"false"`; the last three
restricted to `{"true","false"}`.

The `.cpp` carries about 70 lines of commented-out code - three abandoned
`squared_sum` loops and an abandoned all-pairs `2 * tolerance` sum. None of it
executes and none of it is ported; it is noted here only so a reader comparing
the files does not think something was dropped.

## Preserved source conventions

- **The formula**, `ZhangSimilarityScore.cpp:65-188`: `sum1` and `sum2` are the
  two spectra's **total** intensities - not squared, which is the difference
  from `SpectrumAlignmentScore` - and over every peak pair closer than
  `tolerance`, `sum += sqrt(I1 * I2 * factor)`, giving
  `score = sum / sqrt(sum1 * sum2)` at `:188`.
- **Many-to-many pairing.** There is no alignment: a peak may pair with several
  partners, so a dense region contributes repeatedly and a self-score can exceed
  one (`1.8268` on the upstream fixture).
- **The sliding `j_left` cursor**, `:144-175`. `j_left` is read when a reference
  peak's inner loop starts and written during it, so a write takes effect on the
  *next* reference peak; it advances only when the target peak lies at least the
  tolerance *below* the reference peak, and the loop breaks when the target peak
  is at least the tolerance *above* it. The port reproduces the loop rather than
  an equivalent windowed traversal, including the fact that a matched pair never
  advances the cursor.
- **The window is strict**: `fabs(pos1 - pos2) < tolerance` at `:150`. A pair at
  exactly the tolerance does not count, and it does not break the loop either -
  it advances `j_left` instead, because `pos2 > pos1` is false when the two are
  equal.
- **`I1 * I2` is a `float` multiplication** (`:159`) before it meets the
  `double` factor; see `SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md` for why an `f64`
  product would be exact and therefore wrong.
- **Gaussian factor**: `erfc(mz_difference / (mz_tolerance * 3.0 * sqrt(2.0)))`
  (`:200-201`), reproduced with `libm::erfc` and that grouping.
- **Linear factor**: `(mz_tolerance - mz_difference) / mz_tolerance` (`:204`).
- **The Gaussian wins when both flags are set**, because `:155-157` calls
  `getFactor_(tolerance, diff, use_gaussian_factor)` whenever either is set. The
  sibling `SpectrumAlignmentScore` resolves the clash the other way. Both are
  refused here.
- `is_relative_tolerance` throws `Exception::NotImplemented` at `:62` before any
  work; the parameter is still registered, and so it is here.

## Native differences

- **`is_relative_tolerance` is `Error::Unsupported`**, the mapping of
  `Exception::NotImplemented`. The parameter stays registered so a `Param` tree
  round-trips between C++ and this port; the source's `// TODO remove parameter`
  is recorded here rather than acted on.
- **The Gaussian denominator is evaluated per call.** The source caches it in a
  function-local `static const double` at `:200`, so the first `getFactor_` call
  in a process fixes `mz_tolerance * 3 * sqrt(2)` for the whole process, across
  every instance and every tolerance. Two scorers configured with different
  tolerances then return the same weighting, and which one is right depends on
  call order. Reproducing that would mean a process-global whose value no caller
  can predict; it is recorded as an upstream defect instead and asserted against
  in `zhang_gaussian_scale_follows_the_instance_not_the_first_call`.
- **A zero total intensity on either side yields `Ok(0.0)`.** The source computes
  `0.0 / sqrt(0.0)` and returns NaN. Two empty spectra, one empty spectrum and
  an all-zero spectrum are all this case.
- **Both weighting flags set is an error** rather than a silent Gaussian.
- **Unsorted peaks are refused** with `Error::UnsortedData`. The source does not
  check, although the `j_left` cursor is only correct on sorted input - unsorted
  input yields a silently wrong score, not an error.
- **Negative and non-finite intensities are refused**; a negative product would
  make the radicand negative and the source would return NaN.
- **A non-finite tolerance is refused.**
- **`max_pairs`** bounds the walk, counting candidates as they are examined -
  including those the strict `<` then rejects, because examining them is the cost
  being bounded. Nothing is allocated, so a refusal leaves both inputs untouched.
  The source has no ceiling and its worst case is `|s1| * |s2|`.
- The wave-A `Copy` struct `comparison::ZhangSimilarityScore` remains a separate,
  parameter-free `f64` convenience with a different windowed traversal; it is not
  this port. Unifying them is deferred because the merge would change values
  asserted in `tests/comparison.rs`, which this package may not edit.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| both spectra empty | `Ok(0.0)`; source returns NaN |
| one spectrum empty | `Ok(0.0)`; source returns NaN |
| all intensities zero | `Ok(0.0)`; source returns NaN |
| no peak inside the tolerance | `Ok(0.0)`, as upstream |
| distance exactly the tolerance | not counted; strict `<`, as upstream |
| `is_relative_tolerance` set | `Err(Error::Unsupported)`, the source's `NotImplemented` |
| both weighting flags set | `Err(Error::InvalidValue)`; upstream silently uses the Gaussian |
| `factor` with a zero or non-finite tolerance | `Err(Error::InvalidValue)`; unreachable from `score`, since a pair needs `0 <= d < tolerance` |
| unsorted peaks | `Err(Error::UnsortedData)`; source scores wrongly |
| negative or non-finite intensity | `Err(Error::InvalidValue)` |
| candidates above `max_pairs` | `Err(Error::InvalidValue)`; nothing allocated or mutated |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/ZhangSimilarityScore_test.cpp`
(sha256 `6bd7d6311306518df08ffd6690aa8e7938b9806d5de55fa936fb49f5297264c2`),
six sections, all on `PILISSequenceDB_DFPIANGER_1.dta` normalised with `to_one`.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `ZhangSimilarityScore()` | `zhang_scorer_constructs_with_source_defaults_and_drops` | `name() == "ZhangSimilarityScore"`, `tolerance == 0.2`, three flags false, four parameters | 3 |
| `~ZhangSimilarityScore()` | same | construction then `drop` | 4 |
| `ZhangSimilarityScore(const ZhangSimilarityScore& source)` | `zhang_scorer_copy_and_assignment_carry_name_and_parameters` | `copy.getName()` and `copy.getParameters()` equal the original's, the two upstream `TEST_EQUAL`s | 3 |
| `ZhangSimilarityScore& operator = (...)` | same | assignment over a scorer set to `tolerance = 0.9` restores `0.2` | 3 |
| `double operator () (const PeakSpectrum& spec) const` | `zhang_scorer_reproduces_the_upstream_self_score` | `self_score == 1.8268153570547176`, `TEST_REAL_SIMILAR(score, 1.82682)` | 3 |
| `double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const` | `zhang_scorer_reproduces_the_upstream_pairwise_scores` | `1.82682` again, and truncated to 100 peaks `0.32874865683513527` against `TEST_REAL_SIMILAR(score, 0.328749)` | 3 |

`1.82682` and `0.328749` are tier 3. The full-precision values come from the
same independent Python model of the C++ described in
`SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md`, which reproduces both class-test literals
(tier 4 for the extra digits). The weighting factors are checked against
`math.erfc` and the closed-form linear expression; the resource ceilings, the
degenerate denominators, the strict window boundary and the many-to-many
behaviour are tier 4, derived from the source expressions.
