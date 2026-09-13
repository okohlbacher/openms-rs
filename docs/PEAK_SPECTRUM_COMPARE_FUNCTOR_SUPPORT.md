# PeakSpectrumCompareFunctor support

Port of `src/openms/include/OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h` and
`src/openms/source/COMPARISON/PeakSpectrumCompareFunctor.cpp` at Core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `3e5d1af32c7704e4368398ac5742603ac876405494bb13899216a9f6c1c77fcc`,
`.cpp` sha256 `85dae27b1795883060ea044e49e14d530736ae97b8c485f3f6437b3227b74aa9`).

Rust: [`comparison::PeakSpectrumCompareFunctor`](../src/comparison.rs).
Tests: [`tests/comparison_functors.rs`](../tests/comparison_functors.rs).
Provenance: [`tests/data/comparison_functors_provenance.json`](../tests/data/comparison_functors_provenance.json).

The C++ class is an abstract `DefaultParamHandler` subclass with two pure
virtual call operators and nothing else: the whole header is 6 members, and the
`.cpp` is 4 defaulted or forwarding definitions. It exists so that a caller can
hold a base-class pointer to any spectrum-similarity functor. In Rust that is a
trait, not a struct with a registry; the `#include`s of six derived classes at
the top of the `.cpp` are OpenMS's factory-registration pattern and not a
dependency of the base on its derivatives.

## Port status: partial

**Every member of the header is ported and no type in this crate implements the
trait.** An abstract base exists to be derived from, so a trait with no shipped
implementor is the shape of the port and not the whole of it; the header is
therefore ledgered `partial`, not `complete`. Its sibling
`BinnedSpectrumCompareFunctor.h` is `complete`: all three of *its* source
derivatives ship here.

### The seven derivatives, and where each stands

`PeakSpectrumCompareFunctor.cpp:11-16` `#include`s six derived classes for
factory registration. That is not the full set: `grep "public
PeakSpectrumCompareFunctor" src/openms/include/OpenMS/COMPARISON/*.h` at
`bc9cc12` returns **seven** headers, and `SpectraSTSimilarityScore` is the one
missing from the include list. All seven define the one-spectrum overload as
`operator()(spec, spec)`.

| Derivative | In this crate | Implements the trait | Note |
| --- | --- | --- | --- |
| `SpectrumAlignmentScore` | `comparison::SpectrumAlignmentScore` | no | ported in an earlier wave as a typed `Copy` struct |
| `ZhangSimilarityScore` | `comparison::ZhangSimilarityScore` | no | same |
| `SteinScottImproveScore` | `comparison::SteinScottImproveScore` | no | same |
| `SpectrumPrecursorComparator` | `comparison::SpectrumPrecursorComparator` | no | same |
| `SpectrumCheapDPCorr` | not ported | - | wave B |
| `PeakAlignment` | not ported | - | wave B; the only functor in `COMPARISON/` that registers a `normalized` parameter (`PeakAlignment.cpp:24`) |
| `SpectraSTSimilarityScore` | not ported | - | wave B; not in the `.cpp` include list, and its constructor calls `setName` **without** the `defaultsToParam_()` that the other derivatives call |

The four that exist keep their configuration in typed `Copy` fields rather than
a `DefaultParamHandler`. Giving them one is not a mechanical change: each
registers a real parameter tree upstream, so a handler that carried only the
functor's name would make `handler_mut().set_parameters(...)` silently
ineffective - a second way to configure the functor that the scoring code does
not read. That is the incoherence this wave's fix commit removed from the binned
scores, and it is not worth reintroducing here; wiring the parameters correctly
belongs to the wave that ports those four headers. See the deferrals in the
work-package report.

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class PeakSpectrumCompareFunctor : public DefaultParamHandler` | `pub trait PeakSpectrumCompareFunctor` | a trait, so a functor composes a [`DefaultParamHandler`](../src/param/handler.rs) instead of inheriting one; the trait is object safe, so `&dyn PeakSpectrumCompareFunctor` replaces the base-class pointer |
| (inherited) `DefaultParamHandler` surface | `PeakSpectrumCompareFunctor::handler() -> &DefaultParamHandler`, `handler_mut() -> &mut DefaultParamHandler` | required trait methods; `getName`/`setName`/`getParameters`/`setParameters` are reached through them, and `name()` is a provided shorthand for `handler().name()` |
| `PeakSpectrumCompareFunctor()` | no trait counterpart; each implementor's own constructor | the base constructor only names the handler `"PeakSpectrumCompareFunctor"`, which a derived constructor immediately overwrites with `setName`. The Rust helper `functor_handler(base, name)` performs both steps in that order |
| `PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source)` | `Clone` on the implementor | the C++ copy constructor is `= default`, so it copies the handler and nothing else |
| `~PeakSpectrumCompareFunctor() override` | drop glue | the C++ destructor is `= default`; virtual destruction has no Rust counterpart because a `Box<dyn PeakSpectrumCompareFunctor>` drops through its own vtable |
| `PeakSpectrumCompareFunctor& operator=(const PeakSpectrumCompareFunctor& source)` | assignment of a cloned implementor | the C++ body is the self-assignment guard plus `DefaultParamHandler::operator=`; Rust assignment is unconditional and self-assignment cannot alias |
| `virtual double operator()(const PeakSpectrum& a, const PeakSpectrum& b) const = 0` | `fn score(&self, a: &MSSpectrum, b: &MSSpectrum) -> Result<f64>` | required; returns `Result` so an implementor can report invalid input instead of a wrong or non-finite score |
| `virtual double operator()(const PeakSpectrum& a) const = 0` | `fn self_score(&self, a: &MSSpectrum) -> Result<f64>` | provided, defaulting to `score(a, a)`, which is what all seven source derivatives implement it as |

No member of this header is unported.

## Preserved source conventions

- The class comment's contract - "the value should be greater equal 0" - is
  carried into the trait documentation as a statement about implementors, not as
  a check: the source does not clamp and neither does the trait.
- The two-step naming (`DefaultParamHandler("PeakSpectrumCompareFunctor")` in the
  base, `setName(<derived>)` then `defaultsToParam_()` in the derivative) is
  reproduced by `functor_handler`, so a functor's observable name and its empty
  parameter tree match the source.
- `self_score` delegates to `score(a, a)` rather than to a closed form, because
  all seven source derivatives do exactly that
  (`SpectrumAlignmentScore.cpp:42`, `ZhangSimilarityScore.cpp:47`,
  `SteinScottImproveScore.cpp:50`, `SpectrumPrecursorComparator.cpp:39`,
  `SpectrumCheapDPCorr.cpp:69`, `PeakAlignment.cpp:44`,
  `SpectraSTSimilarityScore.cpp:37`).

## Native differences

- `score` and `self_score` return `Result<f64>`. The source signature cannot
  fail, so a derived functor that meets invalid input either propagates NaN or
  throws from a helper. Making the failure explicit is the crate's convention.
- `handler_mut` exposes the whole mutable handler where C++ exposes individual
  inherited setters. The source's protected `updateMembers_` hook, which
  `DefaultParamHandler` calls after a parameter change, has no trait counterpart:
  a Rust functor that derives typed state from parameters recomputes it in its
  own setter. That is the same decision `DEFAULT_PARAM_HANDLER_SUPPORT.md`
  records for `set_parameters_with`.
- No shipped implementor; see **Port status** above for the seven derivatives
  and why the four that exist in `comparison` do not implement the trait yet.

## Checked boundaries and evidence

- Bounded work: the trait itself performs no work. Each implementor states its
  own ceiling; the three binned functors share [`MAX_COMPARED_BINS`](../src/comparison.rs).
- No panics: `score` returns `Result`; no indexing, division or `sqrt` occurs in
  the trait's provided methods.
- OpenMP: neither this header nor its `.cpp` carries a `#pragma omp`, so there is
  no parallelism gap to record.

### Class-test sections

`src/tests/class_tests/openms/source/PeakSpectrumCompareFunctor_test.cpp`
(sha256 `68ad0e8913ac57c83f7c9301dbe219aa1efdf84cefe731456ff8cebe1008a032`)
declares six sections and every one of them is `NOT_TESTABLE`, with the comment
"pure interface class cannot test this". Zero assertion macros in the file.

| Section | Rust test | Asserted value |
| --- | --- | --- |
| `PeakSpectrumCompareFunctor()` | `peak_spectrum_compare_functor_base_is_a_named_handler_and_a_pair_of_operators` | `functor.name() == "SharedIntensityProduct"` after the base names the handler `"PeakSpectrumCompareFunctor"` and the derivative renames it |
| `PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor& source)` | same | `copy.handler().parameters() == functor.handler().parameters()` |
| `~PeakSpectrumCompareFunctor()` | same | drop glue; the test constructs and drops the functor and the `&dyn` view of it. No observable value, as upstream |
| `PeakSpectrumCompareFunctor& operator=(...)` | same, and `peak_spectrum_compare_functor_handler_is_mutable_as_in_the_source` | after `handler_mut().set_name("renamed")`, `name() == "renamed"` |
| `double operator()(const PeakSpectrum& a, const PeakSpectrum& b) const` | same | `score` of the test functor over `{100,200,300}` and `{200,300}` is `2*3 + 4*5 == 26.0` |
| `double operator()(const PeakSpectrum& a) const` | same | `self_score` is `1 + 4 + 16 == 21.0`, equal to `score(a, a)` |

Evidence tier 4 (independently derived / Rust-only) throughout: the upstream
sections assert nothing, so there is no literal to transcribe. The values above
are computed by hand from the test functor's own definition.
