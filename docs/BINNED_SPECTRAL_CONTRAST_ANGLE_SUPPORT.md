# BinnedSpectralContrastAngle support

Port of `src/openms/include/OpenMS/COMPARISON/BinnedSpectralContrastAngle.h` and
`src/openms/source/COMPARISON/BinnedSpectralContrastAngle.cpp` at Core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `020ad7a63c27bdc711e9c71dcf22c80b2ddf065961a241d4a18c8422041e6a81`,
`.cpp` sha256 `ea854a8917eedb06c6fa44bc0b49bea6673cac53da3f6f260409f1af34d66ac6`).

Rust: [`comparison::BinnedSpectralContrastAngle`](../src/comparison.rs).
Tests: [`tests/comparison_functors.rs`](../tests/comparison_functors.rs).
Provenance: [`tests/data/comparison_functors_provenance.json`](../tests/data/comparison_functors_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class BinnedSpectralContrastAngle : public BinnedSpectrumCompareFunctor` | `pub struct BinnedSpectralContrastAngle` implementing [`BinnedSpectrumCompareFunctor`](BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition instead of inheritance |
| `BinnedSpectralContrastAngle()` | `BinnedSpectralContrastAngle::new() -> Result<Self>` | reproduces `setName("BinnedSpectralContrastAngle")` then `defaultsToParam_()` (`BinnedSpectralContrastAngle.cpp:21-22`) |
| `BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle& source)` | `Clone` | the C++ copy constructor forwards to the base only |
| `~BinnedSpectralContrastAngle() override` | drop glue | C++ destructor is `= default` |
| `BinnedSpectralContrastAngle& operator=(const BinnedSpectralContrastAngle& source)` | assignment of a clone | self-assignment guard plus base assignment |
| `double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override` | `BinnedSpectrumCompareFunctor::score` | see below |
| `double operator()(const BinnedSpectrum& spec) const override` | `BinnedSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` |
| `protected: void updateMembers_() override` | not ported | empty body |
| `protected: double precursor_mass_tolerance_` | not ported | uninitialised, unwritten, unread; see `BINNED_SHARED_PEAK_COUNT_SUPPORT.md` |

`@htmlinclude OpenMS_BinnedSpectralContrastAngle.parameters` documents an empty
parameter section.

## Preserved source conventions

- **It returns a cosine, not an angle.** Despite the class name and the cited
  paper, `BinnedSpectralContrastAngle.cpp:66` computes
  `numerator / sqrt(sum1 * sum2)` and returns it. No `acos` appears anywhere in
  the file. The Rust documentation says so at the type.
- **Grouping.** The denominator is `sqrt(sum1 * sum2)`, not
  `sqrt(sum1) * sqrt(sum2)`, and the source's grouping is kept verbatim because
  the two are not interchangeable. `sum1` is an `f32` value widened to `f64`, so
  `sum1 * sum1` is exact in `f64` - a 24-bit significand squares into 48 bits -
  and `sum1 / sqrt(sum1 * sum1)` is therefore exactly `1.0` for every input.
  Under the other grouping it is not: over 200 000 random `f32`-valued sums,
  `s / (sqrt(s) * sqrt(s))` differs from `1.0` in 93 785 cases, 47% of them. A
  "cleaner" rewrite would have made self-similarity only approximately one.
- **The reduction is `f32`, not `f64`.** `Eigen::SparseVector<float>::dot`
  returns `Scalar`, i.e. `float`; the three `const double` results at
  `BinnedSpectralContrastAngle.cpp:55-57` are widenings of a `float`
  accumulation. The port reproduces that: `sparse_dot` accumulates
  `res += a[i] * b[i]` in `f32` over the shared indices in ascending index order
  - the order Eigen's `SparseDot.h` walk produces - and widens once on return.
  An independent Python reimplementation of both variants over the class-test
  fixture at bin size 1.5, spread 2 and offset 0.4 gives `0.9999813234287416`
  for the `f32` reduction against `0.9999812942748718` for an `f64` one, a
  relative difference of `2.9e-8`. That is small, and it is exactly the kind of
  difference that accumulates silently through a scoring pipeline, so the source
  precision is what the port reproduces.
- **The degenerate guard is the source's own.** `if (sum1 * sum2 == 0.0) return
  0.0;` (`BinnedSpectralContrastAngle.cpp:61-64`), with the source's comment
  naming it a regression fix for `0.0 / 0.0 = NaN`. The port keeps the guard, the
  condition and the value.
- **The score is not clamped.** Bins may be negative in principle, and the source
  neither clamps nor takes an absolute value; a spectrum compared against its own
  negation scores exactly `-1.0`.

## Native differences

- **Incompatible binning is an error, not an assertion.** The source uses
  `OPENMS_PRECONDITION` (`BinnedSpectralContrastAngle.cpp:52`), which is a no-op
  in release builds, so a release-build caller silently scores two different
  binnings against each other. Its two sibling functors throw
  `Exception::IllegalArgument` for the same condition, and their class tests were
  updated with the note "regression: the precondition was a no-op in release
  builds" - this one was not. The port returns `Error::InvalidValue` from all
  three, and the inconsistency is recorded as a C++ issue in the work-package
  report.
- A non-finite dot product - `f32` overflow, reachable at bin intensities around
  `3e38` - is `Error::InvalidValue` rather than an infinity or NaN propagated
  into the score.
- `sqrt` is taken through a helper that refuses a negative or non-finite
  radicand. The radicand here is a product of two sums of squares and so cannot
  be negative; the helper's reachable branch is the non-finite one.
- The comparison is refused above [`MAX_COMPARED_BINS`](../src/comparison.rs)
  combined stored bins. The source has no ceiling.
- The pre-existing free function `comparison::binned_cosine` computes the same
  score in `f64` with a `sqrt(sum1) * sqrt(sum2)` denominator and a clamp to
  `[-1, 1]`. It is retained unchanged for its existing callers; this functor is
  the faithful port of the header, and collapsing the two is a follow-up.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| incompatible unit, size or offset | `Err(Error::InvalidValue)`, both argument orders; source asserts only in debug |
| differing spread | accepted, as upstream |
| empty or all-zero spectrum on either side | `Ok(0.0)`, the source's own guard |
| both empty | `Ok(0.0)` |
| negative bins | permitted and unclamped; `-1.0` against the negation |
| `f32` dot-product overflow | `Err(Error::InvalidValue)` |
| combined stored bins above `MAX_COMPARED_BINS` | `Err(Error::InvalidValue)` before any traversal |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/BinnedSpectralContrastAngle_test.cpp`
(sha256 `6d1557241f32c9504abf76fcfefb72b387820f330fcb540671b1486dac11a03f`),
six sections, over `PILISSequenceDB_DFPIANGER_1.dta`.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `BinnedSpectralContrastAngle()` | `binned_spectral_contrast_angle_constructs_copies_and_assigns` | `name() == "BinnedSpectralContrastAngle"`, parameters empty | 4 |
| `~BinnedSpectralContrastAngle()` | same | construction and drop | 4 |
| `BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle& source)` | same | `copy.handler().parameters() == functor.handler().parameters()` | 3 |
| `BinnedSpectralContrastAngle& operator=(const BinnedSpectralContrastAngle& source)` | same | assignment restores `"BinnedSpectralContrastAngle"` over `"scratch"` | 3 |
| `double operator()(const BinnedSpectrum&, const BinnedSpectrum&) const` | `binned_spectral_contrast_angle_scores_the_upstream_fixture`, `binned_spectral_contrast_angle_is_a_cosine_over_the_bin_vectors`, `binned_functors_reject_incompatible_binning` | `score(bs1, bs2) == 0.999985` for `BinnedSpectrum(s, 1.5, false, 2, 0.4)` against the same spectrum with its last peak dropped; the empty-spectrum regression case scores exactly `0.0` | 3 for the literal, 4 for the guard |
| `double operator()(const BinnedSpectrum&) const` | `binned_spectral_contrast_angle_self_similarity_is_one` | `self_score(bs1) == 1.0` exactly | 4 |

`0.999985` is tier 3, transcribed. Everything else is derived and asserted
exactly rather than to a tolerance: self-similarity is `sum1 / sqrt(sum1 * sum1)`
and therefore exactly `1.0`; orthogonal bin vectors score exactly `0.0`;
`{3, 4}` against `{4, 3}` scores exactly `24/25`, symmetric in both argument
orders; and `{3, 4}` against `{-3, -4}` scores exactly `-1.0`.
