# SpectrumAlignment support

Port of `src/openms/include/OpenMS/COMPARISON/SpectrumAlignment.h` and
`src/openms/source/COMPARISON/SpectrumAlignment.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `e38817397c8f34619b0231887d0b493928f8889163da88ab6eee275e8313e021`,
`.cpp` sha256 `34c39c22df3d5cb7d5d2781a39959023d80084277c40db9431d6ca8a4f3af008`).

Rust: [`comparison::SpectrumAligner`](../src/comparison.rs), which configures
[`comparison::SpectrumAlignment`](../src/comparison.rs) - the banded dynamic
program and the `MatchedIterator` walk - from its parameter tree.
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

This header is the shared primitive of the whole package: its tolerance
semantics and its tie-breaking decide what
[`SpectrumAlignmentScore`](SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md) sees.

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SpectrumAlignment : public DefaultParamHandler` | `pub struct SpectrumAligner` | composition of a `DefaultParamHandler` instead of inheritance. The type name `SpectrumAlignment` was already taken in this module by the wave-A configuration struct this class delegates to; see "Native differences" |
| `SpectrumAlignment()` | `SpectrumAligner::new() -> Result<Self>` | reproduces `SpectrumAlignment.cpp:17-22`: handler named `"SpectrumAlignment"`, `tolerance` `0.3`, `is_relative_tolerance` `"false"` with the `{"true","false"}` restriction, then `defaultsToParam_()`. Fallible only because `DefaultParamHandler::new` is |
| `SpectrumAlignment(const SpectrumAlignment& source)` | `Clone` | the C++ copy constructor is `= default` and copies the base's parameter trees |
| `~SpectrumAlignment() override` | drop glue | C++ destructor is `= default` |
| `SpectrumAlignment& operator=(const SpectrumAlignment& source)` | assignment of a clone | the C++ body is the self-assignment guard plus `DefaultParamHandler::operator=` |
| `template <typename SpectrumType1, typename SpectrumType2> void getSpectrumAlignment(vector<pair<Size, Size>>& alignment, const SpectrumType1& s1, const SpectrumType2& s2) const` | `SpectrumAligner::spectrum_alignment(&self, s1, s2) -> Result<Vec<(usize, usize)>>` | the out-parameter becomes the return value; the source clears it first, so nothing is lost. The two template parameters existed so a `PeakSpectrum` could be aligned against an `MSSpectrum` of another peak type; the port takes two `MSSpectrum` values, which is the only instantiation in the SDK |
| inherited `getParameters` / `setParameters` / `getName` / `setName` | `handler()`, `handler_mut()`, `name()` | `handler_mut().set_parameters(&p)` is `setParameters` |
| - | `SpectrumAligner::tolerance() -> Result<Tolerance>` | native: the parameter pair read back as the typed window the alignment applies |
| - | `SpectrumAligner::max_cells` | native resource ceiling, default [`DEFAULT_ALIGNMENT_CELLS`](../src/comparison.rs) = 5000000 |
| `#define ALIGNMENT_DEBUG` / `#undef ALIGNMENT_DEBUG` at the top of the header | not ported | the macro is defined and immediately undefined, so all five `#ifdef ALIGNMENT_DEBUG` blocks - matrix dumps and alignment printouts to `cerr` - are dead in every build |
| `@htmlinclude OpenMS_SpectrumAlignment.parameters` | the parameter table in the rustdoc | |
| `TODO: improve time complexity, currently O(|s1|*log(|s2|))` | discharged | the ppm path here is a single forward merge over both peak lists, `O(|s1| + |s2|)`, because `MatchedIterator` never restarts its target cursor |

## Preserved source conventions

- **The banded dynamic program**, `SpectrumAlignment.h:79-176`. Gap cost is the
  tolerance, the alignment cost is the m/z distance, and intensity plays no role.
  The row-0 and column-0 initialisation (`i * tolerance`, `j * tolerance`), the
  `left_ptr` band tightening at `:119`, the `off_band` early break at `:112`, the
  three-way minimum with its `<=` tie rules at `:157` and the fallback cost
  `(i - 1 + j - 1) * tolerance` for cells the sparse `std::map` never stored are
  all reproduced. `tests/comparison.rs` carries a map-shaped transcription of
  this loop and cross-checks the compact banded implementation against it over
  200 randomised cases at five tolerances.
- **Traceback starts at the last cell that chose the diagonal**
  (`SpectrumAlignment.h:241`), not at the bottom-right corner, so a trailing run
  of gaps is discarded. Where the C++ reads a `traceback` entry that was never
  written, `std::map::operator[]` default-constructs `(0, 0)` and the
  `while (i >= 1 && j >= 1)` loop ends; the port breaks out at the same point.
- **The ppm path is `MatchedIterator<PpmTrait>`** (`SpectrumAlignment.h:286`),
  with its `float` arithmetic intact: the window is
  `Math::ppmToMass(tolerance, (float)mz)` instantiated at `float`, and the
  distance is `fabs(double - double)` narrowed to `float` on return. The cursor
  into the target never moves backwards across reference peaks, and it stops at
  the first equal adjacent distance - so a duplicate target peak hides a nearer
  one behind it. That quirk is asserted in `tests/comparison.rs`.
- **Asymmetry of the ppm path.** A target peak may serve several reference
  peaks, a reference peak at most one target peak, and the window is derived
  from the *reference* m/z, so `align(a, b)` and `align(b, a)` differ. The
  header's two `@note` lines say exactly this and both are carried into the
  rustdoc.
- Two empty spectra, or one empty spectrum, give an empty alignment upstream and
  here.

## Native differences

- **Unsorted input** is `Error::UnsortedData`, the source's
  `Exception::IllegalArgument` at `:66`. The check is upstream's own, not an
  addition.
- **A negative tolerance** yields an empty alignment upstream - no cell can
  satisfy `diff_align <= tolerance` - and `Error::InvalidValue` here. A negative
  matching window is a configuration mistake, and silently scoring zero hides it.
- **Non-finite coordinates** are rejected. The source would propagate NaN
  through every comparison and return whatever the traceback happened to reach.
- **A ppm coordinate outside the `f32` range** is rejected rather than silently
  becoming an infinity inside `MatchedIterator`'s `float` arithmetic.
- **`max_cells`** bounds the dynamic program. `|s1| + |s2| + 1` is charged before
  the matrix exists and each cell as it is filled, so a refusal allocates at
  most the rows already built and leaves both inputs untouched. The source has no
  ceiling and allocates a `std::map` node per cell.
- **Naming.** The wave-A module already exported a `Copy` configuration struct
  called `SpectrumAlignment` holding a `Tolerance` and `max_cells`; that struct
  *is* this algorithm and this port calls it. The parameterised class is
  therefore `SpectrumAligner`. There is one implementation of the alignment, not
  two.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| both spectra empty | `Ok(vec![])`, as upstream |
| one spectrum empty | `Ok(vec![])`, as upstream |
| no peak inside the tolerance | `Ok(vec![])`; the traceback never leaves `(0, 0)` |
| either spectrum unsorted | `Err(Error::UnsortedData)` |
| negative or non-finite tolerance | `Err(Error::InvalidValue)`; upstream returns an empty alignment |
| non-finite m/z | `Err(Error::InvalidValue)` |
| ppm m/z beyond `f32` | `Err(Error::InvalidValue)`; upstream narrows to an infinity |
| cells above `max_cells` | `Err(Error::InvalidValue)`; nothing is returned and nothing is mutated |
| exact tolerance boundary | `diff <= tolerance` is inclusive in the banded path, `diff <= allowed` in the ppm path, both as upstream |

OpenMP: no `#pragma omp` in this header, its `.cpp` or `MatchedIterator.h`. The
port is serial because the source is.

### Class-test sections

`src/tests/class_tests/openms/source/SpectrumAlignment_test.cpp`
(sha256 `7c31d94741ffd71086ad2dccf4f214c1d8692c1b3033dff6c7759e6f119a0bc1`),
five sections. Fixtures `SpectrumAlignment_in1.dta`, `SpectrumAlignment_in2.dta`
and `PILISSequenceDB_DFPIANGER_1.dta` are retained byte-identical as
`tests/data/comparison_alignment_1.dta`, `tests/data/comparison_alignment_2.dta`
and `tests/data/comparison_dfpianger.dta`.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SpectrumAlignment()` | `spectrum_aligner_constructs_with_source_defaults_and_drops` | `name() == "SpectrumAlignment"`, `tolerance == 0.3`, `is_relative_tolerance == false`, valid strings `["true", "false"]` | 3 |
| `virtual ~SpectrumAlignment()` | same | construction then `drop`; upstream is `delete ptr` | 4 |
| `SpectrumAlignment(const SpectrumAlignment& source)` | `spectrum_aligner_copy_and_assignment_carry_name_and_parameters` | after `tolerance = 0.2`, `copy.name() == first.name()` and `copy.handler().parameters() == first.handler().parameters()` - the upstream `TEST_EQUAL` pair | 3 |
| `SpectrumAlignment& operator=(const SpectrumAlignment& source)` | same | assignment over a default-constructed aligner replaces `0.3` with the assigned tree, `TEST_EQUAL(p, sas2.getParameters())` | 3 |
| `getSpectrumAlignment(...)` (14 assertion macros, so ported rather than mapped) | `spectrum_aligner_reproduces_the_upstream_alignment_golden` | `align(dfpianger, dfpianger).len() == 127`; truncated to 100 peaks, `== 100`; at `tolerance = 1.01`, `[(0,0),(1,1),(3,3),(4,5),(6,6)]`; at 10 ppm, `[(6,6)]`; at 10000 ppm, `[(0,0),(1,1),(2,2),(3,3),(4,5),(5,5),(6,6)]` | 3 |

All five index pairs and all three sizes are transcribed class-test literals
(tier 3). They were also reproduced ahead of the Rust code by an independent
Python transcription of `SpectrumAlignment.h` and of `MatchedIterator`, written
from the C++ rather than from the port; that model agrees with every literal
above, which is what gives confidence the Rust is reading the same algorithm
rather than the same expectations. It remains tier 3/4 evidence: no C++ was
executed, and there is no retained upstream output for this class.
