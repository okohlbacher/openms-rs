# BinnedSumAgreeingIntensities support

Port of `src/openms/include/OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h` and
`src/openms/source/COMPARISON/BinnedSumAgreeingIntensities.cpp` at Core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `5677ec6390ab092c608358df8dfa75830417a3a942d5eb01632384134c7aabd2`,
`.cpp` sha256 `a0f3e2553b3858c2d8dc8cfa5340c2221b208d989c1b7978988ba349bb67fff0`).

Rust: [`comparison::BinnedSumAgreeingIntensities`](../src/comparison.rs).
Tests: [`tests/comparison_functors.rs`](../tests/comparison_functors.rs).
Provenance: [`tests/data/comparison_functors_provenance.json`](../tests/data/comparison_functors_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class BinnedSumAgreeingIntensities : public BinnedSpectrumCompareFunctor` | `pub struct BinnedSumAgreeingIntensities` implementing [`BinnedSpectrumCompareFunctor`](BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition instead of inheritance |
| `BinnedSumAgreeingIntensities()` | `BinnedSumAgreeingIntensities::new() -> Result<Self>` | reproduces `setName("BinnedSumAgreeingIntensities")` then `defaultsToParam_()` (`BinnedSumAgreeingIntensities.cpp:21-22`) |
| `BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities& source)` | `Clone` | the C++ copy constructor forwards to the base only |
| `~BinnedSumAgreeingIntensities() override` | drop glue | C++ destructor is `= default` |
| `BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities& source)` | assignment of a clone | self-assignment guard plus base assignment |
| `double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override` | `BinnedSpectrumCompareFunctor::score` | see below |
| `double operator()(const BinnedSpectrum& spec) const override` | `BinnedSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` |
| `protected: void updateMembers_() override` | not ported | empty body |
| `protected: double precursor_mass_tolerance_` | not ported | uninitialised, unwritten, unread; see `BINNED_SHARED_PEAK_COUNT_SUPPORT.md` |

`@htmlinclude OpenMS_BinnedSumAgreeingIntensities.parameters` documents an empty
parameter section.

## Preserved source conventions

- **The definition of "agreeing" is ported literally.** Per bin,
  `x = (a + b) * 0.5 - |a - b|`, then `y = max(0, x)`, then `sum_nn = sum(y)`,
  and the score is `min(sum_nn / ((sum1 + sum2) / 2.0), 1.0)`
  (`BinnedSumAgreeingIntensities.cpp:58-69`, where the three steps are the
  source's own numbered comments). The class comment states the same rule in
  words: "Bins whose intensity differences are larger than their average
  intensity receive a weight of zero."
- **Over the union, not the intersection.** The source's expression is
  `((bins1 + bins2) * 0.5) - (bins1 - bins2).cwiseAbs()`, and sparse addition
  yields the union of the two stored index sets. The port walks both ordered maps
  as a merge, in ascending index order, visiting union members exactly as Eigen
  does. A bin present in only one spectrum evaluates to `(v + 0)/2 - |v| <= 0`
  and is truncated away, so the union walk and an intersection walk agree
  numerically, but the union is what the source writes and what the port does.
- **The arithmetic is `f32`.** Every operand of the bin expression is a
  `SparseVector<float>` coefficient and the scalar `0.5` is promoted *down* to
  `float` by Eigen's `promote_scalar_arg`, so the coefficients of `s` and the
  reduction `s.coeffs().cwiseMax(0).sum()` are all `float`; likewise the two
  `const double sum1/sum2` are widenings of `SparseVector<float>::sum()`. The
  port accumulates in `f32` throughout and widens only for the final division,
  which is the source's one `double` operation. An independent Python
  reimplementation of both variants over the class-test fixture at bin size 1.5,
  spread 2 and offset 0.0 gives `0.9970718147928658` for the `f32` arithmetic
  against `0.9970718210319204` for an `f64` rewrite, a relative difference of
  `6.3e-9`; at offset 0.4 the difference is `2.8e-8`.
- **The reduction *association* is not the source's, and the fidelity claim is
  bounded accordingly.** Both reductions in this functor go through Eigen's
  *dense* redux, not a sparse walk: `SparseVector::sum()` maps the stored-value
  array to a dense vector and reduces that, and `s.coeffs()` is likewise a dense
  `Map` over the stored-value array, so `s.coeffs().cwiseMax(0).sum()` is a dense
  reduction too. Eigen vectorises those into several packet accumulators combined
  at the end, an association that depends on the target's packet width and on the
  Eigen version, and that is not the sequential one. The port sums the same `f32`
  values in the same bin order sequentially: what it reproduces is the `f32`
  *precision* of the reduction, and agreement with a vectorised C++ build is to
  `f32` reduction rounding rather than bit for bit. Sequential order is the
  choice because it is deterministic, is what a scalar build produces, and is the
  association a parallel reduction would have to reproduce under
  `src/concept/parallel.rs`.
  `tests/comparison_functors.rs::binned_reductions_are_sequential_f32_in_ascending_bin_order`
  pins it on an input where the two associations differ by one `f32` ulp. The
  sibling `BinnedSpectralContrastAngle` is not affected: Eigen's sparse `dot` is
  a scalar merge whose association the port does reproduce.
- **`cwiseMax(0)` semantics.** Eigen's `numext::maxi(a, 0)` keeps `a` unless
  `a < 0`, so `-0.0` survives; the port writes `if value < 0.0 { 0.0 } else { value }`
  for the same reason rather than `f32::max`, whose NaN handling differs.
- **The upper cap.** `min(..., 1.0)` is kept.

## Native differences

- **A zero mean total intensity yields `Ok(0.0)`.** The source divides by
  `(sum1 + sum2) / 2.0` unguarded: two empty spectra or all-zero bins give
  `0.0 / 0.0 = NaN`, and `std::min(NaN, 1.0)` returns NaN because `1.0 < NaN` is
  false. A positive numerator over a zero denominator gives `+inf`, which `min`
  then caps at `1.0`. The port returns a defined `0.0`, matching the guard its
  sibling `BinnedSpectralContrastAngle` already carries for the same situation.
- Incompatible binning is `Error::InvalidValue` rather than
  `Exception::IllegalArgument`. The source's message says "different bin size or
  spread (incompatible binning)", but `BinnedSpectrum::isCompatible` compares
  unit, size and offset and ignores spread; the header's `@throw` text is the
  accurate one and the port follows it.
- `f32` overflow in either total or in the agreeing sum is
  `Error::InvalidValue`, not an infinity carried into the division.
- Refused above [`MAX_COMPARED_BINS`](../src/comparison.rs) combined stored bins.
- **Negative bins are accepted, as upstream.** The consequence is worth stating:
  a wholly negative spectrum scores `0` against *itself*, because
  `(v + v)/2 - 0 = v` is below zero and is truncated, so the class comment's
  "Perfect agreement results in a similarity score of 1.0" holds only for
  nonnegative bins. That is source behaviour, reproduced, not a port choice.
- **`comparison::binned_sum_agreeing_intensities` is now this implementation, not
  a second one.** It used to be a separate `f64` port that walked the
  intersection and *rejected* negative bins, so the module exported two public
  items that answered the same question differently and disagreed both in the
  last bits and on whether a negative bin is an error at all. This functor is the
  faithful port and is authoritative; the function is the same code path for a
  caller that wants neither a parameter handler nor a
  `&dyn BinnedSpectrumCompareFunctor`. The visible consequence of the collapse is
  that a negative bin is no longer an error from either shape - the source
  truncates it away - which
  `tests/comparison.rs::binned_compatibility_signed_cosine_and_symmetric_scores`
  now asserts as `Ok(0.0)` where it previously asserted `Err`.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| incompatible unit, size or offset | `Err(Error::InvalidValue)`, both argument orders |
| differing spread | accepted, as upstream |
| both spectra empty | `Ok(0.0)`; source returns NaN |
| all-zero bins | `Ok(0.0)`; source returns NaN |
| bins that cancel to a zero mean total | `Ok(0.0)`; source returns NaN or `1.0` |
| negative bins | accepted, as upstream; truncation can make self-similarity `0`. Not an error from the functor or from the parameter-free function |
| `f32` overflow of a total or of the agreeing sum | `Err(Error::InvalidValue)` |
| combined stored bins above `MAX_COMPARED_BINS` | `Err(Error::InvalidValue)` before any traversal |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/BinnedSumAgreeingIntensities_test.cpp`
(sha256 `ecaba07a9d49cb7caa330792a3849695590615877ad6316373bc62649586a616`),
six sections, over `PILISSequenceDB_DFPIANGER_1.dta`.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `BinnedSumAgreeingIntensities()` | `binned_sum_agreeing_intensities_constructs_copies_and_assigns` | `name() == "BinnedSumAgreeingIntensities"`, parameters empty | 4 |
| `~BinnedSumAgreeingIntensities()` | same | construction and drop | 4 |
| `BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities& source)` | same | `copy.handler().parameters() == functor.handler().parameters()` | 3 |
| `BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities& source)` | same | assignment restores `"BinnedSumAgreeingIntensities"` over `"scratch"` | 3 |
| `double operator()(const BinnedSpectrum&, const BinnedSpectrum&) const` | `binned_sum_agreeing_intensities_scores_the_upstream_fixture`, `binned_sum_agreeing_intensities_discards_bins_that_disagree`, `binned_functors_reject_incompatible_binning` | `score(bs1, bs2) == 0.99707` for `BinnedSpectrum(s, 1.5, false, 2, 0.0)` against the same spectrum with its last peak dropped; `score(bs1, bs1) == 1.0` exactly; a bin size of 2.0 yields `Err` | 3 for the literal, 4 for the rest |
| `double operator()(const BinnedSpectrum&) const` | `binned_sum_agreeing_intensities_self_similarity_is_one` | `self_score(bs1) == 1.0` exactly at offset `0.4` | 4 |

`0.99707` is tier 3, transcribed. The rest are derived from the formula and
asserted exactly: self-similarity is `1.0` because `(v + v) * 0.5` is exact in
`f32` and `|v - v| = 0`, so the numerator accumulates the same values in the same
order as `sum1` and the denominator is `(sum1 + sum1) / 2 = sum1`;
`4` against `2` scores `(3 - 2) / 3 = 1/3` in both argument orders; `4` against
`1` scores `0` because `|4 - 1| = 3 > 2.5`; a bin present in only one spectrum
scores `0`; and an all-negative spectrum scores `0` against itself.
