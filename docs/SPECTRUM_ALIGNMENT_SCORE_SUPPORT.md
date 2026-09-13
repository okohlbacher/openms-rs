# SpectrumAlignmentScore support

Port of `src/openms/include/OpenMS/COMPARISON/SpectrumAlignmentScore.h` and
`src/openms/source/COMPARISON/SpectrumAlignmentScore.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `25ebe09fe97a7503923968842ce5939fbcb1ab9894165d58255bd11370632404`,
`.cpp` sha256 `2e08a1e8650f3b071068ab0094682ad5b01f2ac0dddcad8cc8c98d97cf8ea111`).

Rust: [`comparison::SpectrumAlignmentScorer`](../src/comparison.rs), implementing
[`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SpectrumAlignmentScore : public PeakSpectrumCompareFunctor` | `pub struct SpectrumAlignmentScorer` implementing `PeakSpectrumCompareFunctor` | composition of a `DefaultParamHandler` instead of inheritance. This is the first implementor of that trait in the crate |
| `SpectrumAlignmentScore()` | `SpectrumAlignmentScorer::new() -> Result<Self>` | reproduces `SpectrumAlignmentScore.cpp:15-27`: base handler `"PeakSpectrumCompareFunctor"`, `setName("SpectrumAlignmentScore")` at `:18`, four defaults, `defaultsToParam_()` |
| `SpectrumAlignmentScore(const SpectrumAlignmentScore& source)` | `Clone` | C++ copy constructor is `= default` |
| `~SpectrumAlignmentScore() override` | drop glue | C++ destructor is `= default` |
| `SpectrumAlignmentScore& operator=(const SpectrumAlignmentScore& source)` | assignment of a clone | C++ body is the self-assignment guard plus the base assignment |
| `double operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const override` | `PeakSpectrumCompareFunctor::score` | `Result<f64>`; see below |
| `double operator()(const PeakSpectrum& spec) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override at `:42-45` is `return operator()(spec, spec)` |
| inherited `getParameters` / `setParameters` / `getName` | `handler()`, `handler_mut()`, `name()` | |
| - | `SpectrumAlignmentScorer::max_cells` | native resource ceiling forwarded to the alignment |
| `@htmlinclude OpenMS_SpectrumAlignmentScore.parameters` | the parameter table in the rustdoc | |

Parameters: `tolerance` `0.3`, `is_relative_tolerance` `"false"`,
`use_linear_factor` `"false"`, `use_gaussian_factor` `"false"`; the last three
restricted to `{"true","false"}`.

## Preserved source conventions

- **The formula**, `SpectrumAlignmentScore.cpp:68-99`:
  `sum1 = sum pow(I1, 2)`, `sum2 = sum pow(I2, 2)` over all peaks, and over the
  aligned pairs `sum += sqrt(I1 * I2 * factor)`, giving
  `score = sum / sqrt(sum1 * sum2)`.
- **The denominator grouping is one `sqrt` of a product** (`:99`), not the
  product of two `sqrt`s. The two differ in the last bits and the source's form
  is reproduced.
- **`I1 * I2` is a `float` multiplication.** `Peak1D::getIntensity()` returns
  `float`, so `s1[..].getIntensity() * s2[..].getIntensity()` rounds to `float`
  before it meets the `double` factor at `:96`. An `f64` product of two `f32`
  values would be *exact* and therefore wrong: on the upstream `DFPIANGER`
  fixture the two differ by 3.5e-9 relative (1.484501010820065 against
  1.4845010143546342). The port reproduces the `float` step.
- **`pow(getIntensity(), 2)` is the `double` overload** (`:71`), so each
  intensity is widened first and the square is exact; `f64::from(i).powi(2)` is
  the same single multiplication.
- **Gaussian factor**: `epsilon = mz_difference / (3.0 * mz_tolerance * sqrt(2))`
  then `std::erfc(epsilon)` (`:91-92`), reproduced with `libm::erfc` and that
  exact grouping.
- **Linear factor**: `(mz_tolerance - mz_difference) / mz_tolerance` (`:87`).
- **Linear wins when both flags are set**, because `:85-92` is an `if` / `else
  if`. The sibling `ZhangSimilarityScore` resolves the same clash the other way.
  Both are refused here; see below.
- **Two windows, not one.** The alignment is delegated to a fresh
  `SpectrumAlignment` carrying only `tolerance` and `is_relative_tolerance`
  (`:56-60`), which under ppm matches through `MatchedIterator`'s `float`
  `Math::ppmToMass(tol, (float)mz)`. The score then re-derives its own window at
  `:81` as `tolerance * s1[ap.first].getMZ() * 1e-6`, in `double` and with a
  different grouping. The two disagree, so a matched pair can have
  `mz_difference > mz_tolerance` and a negative linear factor. Both expressions
  are reproduced exactly, including the disagreement.
- An alignment with no pairs scores `0.0 / sqrt(sum1 * sum2)`, an unremarkable
  `0.0`, upstream and here.

## Native differences

- **A zero squared-intensity sum on either side yields `Ok(0.0)`.** The source
  computes `0.0 / sqrt(0.0)` and returns NaN. This is the case for two empty
  spectra, for one empty spectrum and for a spectrum whose peaks are all zero.
  The module's binned scorers already make the same choice, and a NaN similarity
  silently poisons every consumer.
- **Both weighting flags set is an error.** Upstream it is
  `OPENMS_PRECONDITION` (`:54`), compiled out of a release build, after which the
  linear factor silently wins.
- **A negative radicand is an error.** `sqrt` of a negative weighted product
  returns NaN upstream; it is reachable through the ppm window mismatch above
  and is reported as `Error::InvalidValue`.
- **A zero `mz_tolerance` under a weighting flag is an error** rather than a
  `0/0` NaN. With `tolerance = 0` only exactly-coincident peaks align, so the
  source evaluates `(0 - 0) / 0`.
- **Unsorted peaks, non-finite values and negative intensities are refused.**
  The source checks sortedness inside the alignment only; negative intensities
  would make the radicand negative.
- **`max_cells`** bounds the delegated alignment; the source has no ceiling.
- The wave-A `Copy` struct `comparison::SpectrumAlignmentScore` remains a
  separate, parameter-free `f64` convenience: it accumulates the intensity
  product in `f64` and divides by `sqrt(sum1) * sqrt(sum2)`. It is *not* this
  port and its own documentation says so. Collapsing the two onto one
  implementation is the obvious follow-up and is deferred here only because the
  merge would change values asserted in `tests/comparison.rs`, which this
  package may not edit.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| both spectra empty | `Ok(0.0)`; source returns NaN |
| one spectrum empty | `Ok(0.0)`; source returns NaN |
| all intensities zero | `Ok(0.0)`; source returns NaN |
| no peak inside the tolerance | `Ok(0.0)`, as upstream |
| `use_linear_factor` and `use_gaussian_factor` both set | `Err(Error::InvalidValue)`; upstream precondition is a release-build no-op |
| ppm pair weighted outside its own window | `Err(Error::InvalidValue)` from the negative radicand; source returns NaN |
| zero tolerance with a weighting flag | `Err(Error::InvalidValue)`; source divides by zero |
| unsorted peaks | `Err(Error::UnsortedData)` |
| negative or non-finite intensity | `Err(Error::InvalidValue)` |
| alignment beyond `max_cells` | `Err(Error::InvalidValue)` before any score is accumulated |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/SpectrumAlignmentScore_test.cpp`
(sha256 `ca3f64a1af4cce8944224e9b22a2d6cdcf42aa85a1fa018ee98dee4534ab560d`),
six sections, all using `PILISSequenceDB_DFPIANGER_1.dta` normalised with
`Normalizer`'s `to_one`, retained here as `tests/data/comparison_dfpianger.dta`.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SpectrumAlignmentScore()` | `spectrum_alignment_scorer_constructs_with_source_defaults_and_drops` | `name() == "SpectrumAlignmentScore"`, `tolerance == 0.3`, three flags false, four parameters | 3 |
| `virtual ~SpectrumAlignmentScore()` | same | construction then `drop` | 4 |
| `SpectrumAlignmentScore(const SpectrumAlignmentScore&)` | `spectrum_alignment_scorer_copy_and_assignment_carry_name_and_parameters` | `copy.name() == first.name()` and equal parameter trees after `tolerance = 0.2` | 3 |
| `SpectrumAlignmentScore& operator=(...)` | same | assignment replaces the default `0.3` tree | 3 |
| `double operator()(const PeakSpectrum&, const PeakSpectrum&) const` | `spectrum_alignment_scorer_reproduces_the_upstream_pairwise_scores` | `score(s, s) == 1.484501010820065` (`TEST_REAL_SIMILAR(score, 1.48268)` at `TOLERANCE_ABSOLUTE(0.01)`); truncated to 100 peaks, `3.8247227487373805` against `TEST_REAL_SIMILAR(score, 3.82472)` | 3 |
| `double operator()(const PeakSpectrum&) const` | `spectrum_alignment_scorer_self_score_is_the_pairwise_score` | `self_score == score(s, s) == 1.484501010820065`, and `> 1.0` | 3 |

The two class-test literals `1.48268` and `3.82472` are tier 3. The
full-precision values asserted alongside them come from an independent Python
model of `SpectrumAlignment.h` and `SpectrumAlignmentScore.cpp` written from the
C++ before the Rust existed; that model reproduces both class-test literals
inside their `TOLERANCE_ABSOLUTE(0.01)` bands and pins the remaining digits.
That is tier 4 for the extra digits, and it is what makes the `float` intensity
product observable at all: the same model in `f64` gives `1.4845010143546342`,
which is still `TEST_REAL_SIMILAR` to `1.48268` and would have hidden the
divergence. Weighting factors, resource ceilings, degenerate denominators and
the ppm window mismatch are tier 4, derived from the expressions.
