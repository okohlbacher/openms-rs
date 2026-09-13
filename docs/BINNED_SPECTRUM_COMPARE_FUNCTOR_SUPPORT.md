# BinnedSpectrumCompareFunctor support

Port of `src/openms/include/OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h` and
`src/openms/source/COMPARISON/BinnedSpectrumCompareFunctor.cpp` at Core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `9842d1cb14dc7d168bd411e06a4b480d769979d50181b8292bd1daba5c62aa6a`,
`.cpp` sha256 `f0b555a2d76bd65fd79de4ed9c0f268e73a93055b87ca3dd028f6c066531a187`).

Rust: [`comparison::BinnedSpectrumCompareFunctor`](../src/comparison.rs).
Tests: [`tests/comparison_functors.rs`](../tests/comparison_functors.rs).
Provenance: [`tests/data/comparison_functors_provenance.json`](../tests/data/comparison_functors_provenance.json).

The second abstract base of the spectrum-similarity hierarchy, identical in
shape to [`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md)
but over [`BinnedSpectrum`](COMPARISON_SUPPORT.md) rather than `PeakSpectrum`.
Its `.cpp` includes its three derived classes for factory registration; that is
not a dependency of the base on them.

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class BinnedSpectrumCompareFunctor : public DefaultParamHandler` | `pub trait BinnedSpectrumCompareFunctor` | a trait; object safe, so `&dyn BinnedSpectrumCompareFunctor` replaces the base-class pointer |
| (inherited) `DefaultParamHandler` surface | `handler() -> &DefaultParamHandler`, `handler_mut() -> &mut DefaultParamHandler`, `name() -> &str` | as for the peak base |
| `BinnedSpectrumCompareFunctor()` | each implementor's constructor, through the private `functor_handler("BinnedSpectrumCompareFunctor", <derived name>)` | the base names the handler after itself and the derivative overwrites that name; both steps are reproduced in order |
| `BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source)` | `Clone` on the implementor | C++ copy constructor is `= default` |
| `~BinnedSpectrumCompareFunctor() override` | drop glue | C++ destructor is `= default` |
| `BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor& source)` | assignment of a cloned implementor | the C++ body is the self-assignment guard plus `DefaultParamHandler::operator=` |
| `virtual double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const = 0` | `fn score(&self, spec1: &BinnedSpectrum, spec2: &BinnedSpectrum) -> Result<f64>` | required; `Result` instead of an unconditional `double` |
| `virtual double operator()(const BinnedSpectrum& spec) const = 0` | `fn self_score(&self, spec: &BinnedSpectrum) -> Result<f64>` | provided, defaulting to `score(spec, spec)`, which is what all three source derivatives implement it as |
| `private:` (empty access block) | nothing | the header opens an empty `private:` section before `public:`; no member follows it |

No member of this header is unported.

## Preserved source conventions

- The two-step naming and `defaultsToParam_()` on an empty defaults tree, so
  every derived functor reports its own name and an empty `getParameters()`.
- `self_score` delegates to `score(spec, spec)`, matching
  `BinnedSharedPeakCount.cpp:41`, `BinnedSpectralContrastAngle.cpp:41` and
  `BinnedSumAgreeingIntensities.cpp:41`.
- Binning compatibility is the derivatives' business, not the base's, exactly as
  upstream: the base declares no compatibility check and the trait requires none.

## Native differences

- `Result<f64>` return type, as for the peak base.
- **The class comment is wrong about this hierarchy and the port says so.** It
  states that "Functors normalized in the range [0,1] are identifiable at the set
  `normalized` parameter of the ParameterHandler". No binned functor in the
  pinned revision registers a `normalized` parameter - the only functor that does
  is `PeakAlignment` (`PeakAlignment.cpp:25`), which derives from the *other*
  base. All three binned functors are in fact normalised for nonnegative bins,
  and each states that in its own documentation instead.
- The header declares no nested exception type. The class test still has two
  sections for a `BinnedSpectrumCompareFunctor::IncompatibleBinning` nested class
  that no longer exists; see below.

## Checked boundaries and evidence

- Bounded work: none in the trait. Implementors share [`MAX_COMPARED_BINS`](../src/comparison.rs),
  checked before any traversal, and no scorer allocates, so a refusal leaves both
  inputs untouched.
- No panics: `score` returns `Result`; the provided methods contain no indexing,
  division or `sqrt`.
- OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/BinnedSpectrumCompareFunctor_test.cpp`
(sha256 `dc95571fcba4c5532cc52f962b83833bed06636805963f25a35bef542b98b9e1`)
declares eight sections, all `NOT_TESTABLE`, with the comment "interface class
is not testable". Zero assertion macros in the file.

| Section | Rust test | Asserted value |
| --- | --- | --- |
| `BinnedSpectrumCompareFunctor()` | `binned_spectrum_compare_functor_base_is_a_named_handler_and_a_pair_of_operators` | each of the three derivatives reports its own name, e.g. `"BinnedSharedPeakCount"`, with `handler().parameters().is_empty()` |
| `~BinnedSpectrumCompareFunctor()` | same | drop glue; the three functors and their `&dyn` views are constructed and dropped. No observable value, as upstream |
| `BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source)` | `binned_shared_peak_count_constructs_copies_and_assigns` and its two siblings | `copy.handler().parameters() == functor.handler().parameters()` |
| `BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor& source)` | same three | after assigning over a handler renamed to `"scratch"`, `assigned.name() == "BinnedSharedPeakCount"` |
| `virtual double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const = 0` | `binned_spectrum_compare_functor_base_is_a_named_handler_and_a_pair_of_operators` | called through `&dyn`, all three score the DFPIANGER fixture against itself as exactly `1.0` |
| `virtual double operator()(const BinnedSpectrum& spec) const = 0` | same | `functor.self_score(&spectrum) == functor.score(&spectrum, &spectrum)` for all three |
| `[BinnedSpectrumCompareFunctor::IncompatibleBinning] IncompatibleBinning(const char*, int, const char*, const char* message="compared spectra have different settings in binsize and/or binspread")` | `binned_functors_reject_incompatible_binning` | a different bin size and a different offset both yield `Err(Error::InvalidValue(_))` from all three functors, in both argument orders |
| `[BinnedSpectrumCompareFunctor::IncompatibleBinning] virtual ~IncompatibleBinning()` | same | the same error path; the nested type does not exist in the pinned header |

The last two sections name a nested exception class that the pinned header does
not declare. `START_SECTION` stringifies its argument rather than compiling it,
so these sections still pass upstream while testing nothing and documenting an
API that was removed; the incompatible-binning contract now lives in
`Exception::IllegalArgument` thrown from two of the three derivatives. This is
recorded as a C++ issue in the work-package report.

Evidence tier 4 for the trait-shape assertions (no upstream literal exists) and
tier 4 for the self-similarity value `1.0`, which is derived rather than
transcribed: shared bins over the larger bin count is `n/n`; the contrast angle
is `sum1 / sqrt(sum1 * sum1)`, exact because the square of an `f32`-valued `f64`
is exact in `f64`; and agreeing intensity against itself is `(v + v)/2 - 0 = v`
per bin, so the numerator is the bin sum and the denominator is the same value.
