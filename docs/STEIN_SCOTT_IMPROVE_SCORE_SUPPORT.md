# SteinScottImproveScore support

Port of `src/openms/include/OpenMS/COMPARISON/SteinScottImproveScore.h` and
`src/openms/source/COMPARISON/SteinScottImproveScore.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `b1a807a03f0d6c4242b5dd3d3e42b4b764146dda83c3b573a2283806ec06c52c`,
`.cpp` sha256 `a6f2b273987ef23d86a12f70cfb2c45cf3ad7c616b57a1d623af07684bfeb84d`).

Rust: [`comparison::SteinScottImproveScorer`](../src/comparison.rs), implementing
[`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SteinScottImproveScore : public PeakSpectrumCompareFunctor` | `pub struct SteinScottImproveScorer` implementing `PeakSpectrumCompareFunctor` | composition instead of inheritance |
| `SteinScottImproveScore()` | `SteinScottImproveScorer::new() -> Result<Self>` | reproduces `SteinScottImproveScore.cpp:17-24`: base handler, `setName("SteinScottImproveScore")` at `:20`, two defaults, `defaultsToParam_()` |
| `SteinScottImproveScore(const SteinScottImproveScore& source)` | `Clone` | C++ copy constructor is `= default` |
| `~SteinScottImproveScore() override` | drop glue | C++ destructor is `= default` |
| `SteinScottImproveScore& operator=(const SteinScottImproveScore& source)` | assignment of a clone | |
| `double operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const override` | `PeakSpectrumCompareFunctor::score` | `Result<f64>`. The `@brief` "Similarity pairwise score" and its body text are carried into the rustdoc |
| `double operator()(const PeakSpectrum& spec) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override at `:51-54` is `return operator()(spec, spec)`. Its `@param[in] spec` and `@see SteinScottImproveScore()` are folded into the trait method's prose |
| inherited `getParameters` / `setParameters` / `getName` | `handler()`, `handler_mut()`, `name()` | |
| - | `SteinScottImproveScorer::max_pairs` | native resource ceiling, default [`DEFAULT_SCORED_PAIRS`](../src/comparison.rs) = 5000000 |
| `@htmlinclude OpenMS_SteinScottImproveScore.parameters` | the parameter table in the rustdoc | |

Parameters: `tolerance` `0.2` ("defines the absolute error of the mass
spectrometer") and `threshold` `0.2` ("if the calculated score is smaller than
the threshold, a zero is given back"). Neither carries a restriction upstream,
and neither does here - a difference from the other three headers in this
package, and the reason `parameters.valid_strings("tolerance")` is an error.

## Preserved source conventions

- **The formula**, `SteinScottImproveScore.cpp:65-118`. With
  `epsilon = tolerance` and `constant = epsilon / 10000`:
  `sum1`/`sum2` are squared-intensity sums, `sum3`/`sum4` total intensities,
  `z = constant * (sum3 * sum4)`, `sum` adds `I1 * I2` over every peak pair with
  `|d| <= 2 * epsilon`, and `score = (sum - z) / sqrt(sum1 * sum2)`.
- **The window is `2 * epsilon` and its boundary is inclusive** (`:95`), unlike
  `ZhangSimilarityScore`'s strict `< tolerance`. The class comment says "within
  the given mass-to-charge range" without naming the factor of two; the code is
  authoritative and the rustdoc states it.
- **The grouping of `z`** is `constant * (sum3 * sum4)` (`:88`) - the two totals
  multiplied first. Rearranging it changes the last bits.
- **`sum1 += temp * temp` with `double temp = it1.getIntensity()`** (`:76-79`):
  the `float` intensity is widened *before* squaring, so the square is exact.
- **`sum += s1[i].getIntensity() * s2[j].getIntensity()`** (`:97`) is a `float`
  multiplication whose result is then added to a `double`. This is the one place
  in this class where single precision is observable; see
  `SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md`.
- **The sliding `j_left` cursor** (`:89-111`) is the same loop shape as
  `ZhangSimilarityScore`'s and is reproduced identically, including the fact
  that a matched pair never advances the cursor and that `j_left` takes effect
  only from the next reference peak.
- **The threshold comparison narrows to `float`** (`:115`):
  `score < (float)param_.getValue("threshold")`. The `float` is widened again for
  the comparison, so the default `0.2` is compared against `0.20000000298023224`,
  not `0.2`. Reproduced with `to_f32()` followed by `f64::from`.
- Two spectra with no peak inside the window give `sum == 0` and therefore
  `-z / sqrt(sum1 * sum2)`, a small negative number that the default threshold
  reports as `0.0` - upstream and here. A negative threshold exposes the raw
  value, which the port asserts.

## Native differences

- **A zero squared-intensity sum on either side yields `Ok(0.0)`.** The source
  computes `(0 - 0) / sqrt(0)`, gets NaN, finds `NaN < threshold` false and
  returns the NaN. Two empty spectra, one empty spectrum and an all-zero
  spectrum are all this case, and the module's binned scorers already choose a
  defined zero for their own degenerate denominators.
- **Unsorted peaks are refused** with `Error::UnsortedData`. The source does not
  check, although the `j_left` cursor is only correct on sorted input.
- **Negative intensities are refused.** The formula itself tolerates them -
  nothing takes a square root of an intensity here - but the module's scoring
  validation is shared and consistent, and a negative intensity is a data defect
  rather than a modelling choice. This is the one place in this port that is
  stricter than the source without a numerical reason, and it is stated here so
  it is a decision rather than an accident.
- **Non-finite intensities, tolerances and thresholds are refused.**
- **`max_pairs`** bounds the walk, counting candidates as they are examined.
  Nothing is allocated, so a refusal leaves both inputs untouched. The source
  has no ceiling and its worst case is `|s1| * |s2|`.
- The wave-A `Copy` struct `comparison::SteinScottImproveScore` remains a
  separate, parameter-free `f64` convenience with a different traversal and a
  different grouping of `z`; it is not this port. Unifying them is deferred
  because the merge would change values asserted in `tests/comparison.rs`, which
  this package may not edit.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| both spectra empty | `Ok(0.0)`; source returns NaN |
| one spectrum empty | `Ok(0.0)`; source returns NaN |
| all intensities zero | `Ok(0.0)`; source returns NaN |
| no peak inside `2 * tolerance` | raw score `-z / sqrt(sum1 * sum2)`, reported as `0.0` by the default threshold, as upstream; a negative threshold returns the raw value |
| distance exactly `2 * tolerance` | counted; the comparison is `<=`, as upstream |
| default threshold | compared as `0.20000000298023224`, the `float` narrowing of `0.2` |
| unsorted peaks | `Err(Error::UnsortedData)`; source scores wrongly |
| negative or non-finite intensity | `Err(Error::InvalidValue)`; stricter than the source, by choice |
| non-finite tolerance or threshold | `Err(Error::InvalidValue)` |
| candidates above `max_pairs` | `Err(Error::InvalidValue)`; nothing allocated or mutated |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/SteinScottImproveScore_test.cpp`
(sha256 `e385fde342aaedf136d08e59c2a97249487c6d1a1acc625177b7e438b72b1a62`),
six sections. Its fixture is built in the test itself: five peaks at m/z 500,
600, 700, 800 and 900 whose intensity equals their m/z. No file is needed.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SteinScottImproveScore()` | `stein_scott_scorer_constructs_with_source_defaults_and_drops` | `name() == "SteinScottImproveScore"`, `tolerance == 0.2`, `threshold == 0.2`, two parameters, no valid-string restriction | 3 |
| `virtual ~SteinScottImproveScore()` | same | construction then `drop` | 4 |
| `SteinScottImproveScore(const SteinScottImproveScore& source)` | `stein_scott_scorer_copy_and_assignment_carry_name_and_parameters` | `copy.getName()` and `copy.getParameters()` equal the original's | 3 |
| `SteinScottImproveScore& operator = (...)` | same | assignment over a scorer set to `threshold = 0.9` restores the default tree | 3 |
| `double operator () (const PeakSpectrum& spec) const` | `stein_scott_scorer_reproduces_the_upstream_self_score` | `self_score == (2550000 - 0.2 / 10000 * 3500^2) / 2550000 == 0.9999039215686274`, and `> 0.99` as upstream requires before it rounds to 1 | 4, with 3 |
| `double operator () (const PeakSpectrum& spec1, const PeakSpectrum& spec2) const` | `stein_scott_scorer_reproduces_the_upstream_pairwise_score` | the same value on two separately built copies of the fixture | 4, with 3 |

The class test only asserts `if (score > 0.99) score = 1; TEST_REAL_SIMILAR(score,
1)`, which is weak: it passes for anything above `0.99`. The Rust test asserts
the exact value and **derives it** instead of transcribing it, which is stronger
than the oracle. The peaks are 100 Th apart and the window is `+-0.4`, so the
only pairs are each peak with itself; every `float` product `500*500` through
`900*900` is exactly representable, so `sum` is exactly `2550000`, `sum1` and
`sum2` are exactly `2550000`, `sqrt(2550000 * 2550000)` is exactly `2550000`,
and `z = (0.2 / 10000) * 3500^2`. The test asserts that closed form and the
resulting literal. That derivation is tier 4; the `> 0.99` band it lands in is
the tier-3 class-test evidence.
