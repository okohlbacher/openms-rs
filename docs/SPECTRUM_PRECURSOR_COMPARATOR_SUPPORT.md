# SpectrumPrecursorComparator support

Port of `src/openms/include/OpenMS/COMPARISON/SpectrumPrecursorComparator.h` and
`src/openms/source/COMPARISON/SpectrumPrecursorComparator.cpp` at Core SDK
revision `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `bae8cbf781b217b6405853f4d82ec7d4ebcfcb633b4250c1987bd19c42514e48`,
`.cpp` sha256 `eec1ce3b572a00631f3198c1c6717928c136cc93a2c4d1aff3621508d32b2527`).

Rust: [`comparison::SpectrumPrecursorComparator`](../src/comparison.rs).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

Every public member of the header appears here.

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SpectrumPrecursorComparator : public PeakSpectrumCompareFunctor` | `pub struct SpectrumPrecursorComparator` implementing [`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition of a `DefaultParamHandler` instead of inheritance |
| `SpectrumPrecursorComparator()` | `SpectrumPrecursorComparator::new() -> Result<Self>`, and `Default` | `new` is fallible only through the handler's resource limits, which a two-word name and one default entry cannot reach; `Default` calls it and `expect`s that. Reproduces `setName("SpectrumPrecursorComparator")`, `defaults_.setValue("window", 2, ...)`, `defaultsToParam_()` (`SpectrumPrecursorComparator.cpp:21-23`) |
| `SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source)` | `Clone` | the C++ copy constructor is `= default` |
| `~SpectrumPrecursorComparator() override` | drop glue | C++ destructor is `= default` |
| `SpectrumPrecursorComparator& operator=(const SpectrumPrecursorComparator& source)` | `Clone::clone_from`, or assignment of a clone | the C++ body is the self-assignment guard plus the base assignment |
| `double operator()(const PeakSpectrum& a, const PeakSpectrum& b) const override` | `PeakSpectrumCompareFunctor::score` | returns `Result<f64>`; see below |
| `double operator()(const PeakSpectrum& a) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` (`SpectrumPrecursorComparator.cpp:39-42`) |
| the registered `window` parameter | `window() -> Result<f64>`, read through `handler()` / `handler_mut()` | read-only accessor; the source has no setter either, and the parameter is the only way in |

`@htmlinclude OpenMS_SpectrumPrecursorComparator.parameters` documents the one
`window` entry, registered as the **integer** `2` with the description
"Allowed deviation between precursor peaks." and cast to `double` on every call.

## Preserved source conventions

- **Formula.** `if (fabs(mz1 - mz2) > window) return 0; return window - fabs(mz1 - mz2);`
  (`SpectrumPrecursorComparator.cpp:57-62`). That is `max(0, window - |Δ|)`
  without the redundant second subtraction, and the port computes it that way;
  the two agree bit for bit, because the clamp only replaces a value the first
  branch would have returned as an exact `0`.
- **Only the first precursor is read**, and a spectrum with no precursor
  contributes m/z `0.0` (`SpectrumPrecursorComparator.cpp:45-55`). Two spectra
  that both lack a precursor therefore score the full `window`, which is not a
  meaningful similarity but is the source's behaviour and is tested.
- The parameter is re-read from the handler on every call, as
  `param_.getValue("window")` is.

## Native differences

- **The type changed shape in this wave.** An earlier wave shipped
  `SpectrumPrecursorComparator` as a `Copy` struct with a public `window: f64`
  field and an inherent `score`. It is now the header's functor: a
  `DefaultParamHandler` carrying the registered `window` parameter, implementing
  `PeakSpectrumCompareFunctor`. `Default` and the `score` entry point are
  unchanged for callers; the public field is gone, replaced by `window()` over
  the parameter, because a field and a parameter that can disagree is exactly
  what made `set_parameters` ineffective on the earlier shape.
- A `window` that is negative or not finite is `Err(Error::InvalidValue)`. The
  source reads the parameter unchecked and would return a negative similarity
  for every pair.
- `score` returns `Result<f64>`; an invalid spectrum (a non-finite coordinate, a
  malformed precursor) is an error rather than arithmetic on a NaN.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| neither spectrum has a precursor | `Ok(window)`; both m/z read as `0.0`, as upstream |
| one spectrum has no precursor | `Ok(max(0, window - |other m/z|))`, as upstream |
| several precursors | only the first is read, as upstream |
| distance beyond the window | `Ok(0.0)` exactly |
| negative or non-finite `window` | `Err(Error::InvalidValue)`; source returns negative scores |
| invalid spectrum | `Err`; source has no validation |
| peak count | irrelevant: no peak is read, so there is nothing to bound |

OpenMP: no `#pragma omp` in this header or its `.cpp`; the port is serial, as the
source is.

### Class-test sections

`src/tests/class_tests/openms/source/SpectrumPrecursorComparator_test.cpp`
(sha256 `15c1c5e6054344e0d2f54c6fe06f5898a3e71bc795f72d3eccf153b029f30a03`),
**six sections, all ported**. The fixtures are `Transformers_tests.dta` and
`Transformers_tests_2.dta`, retained byte-identical as
[`tests/data/Transformers_tests.dta`](../tests/data/Transformers_tests.dta) and
[`tests/data/comparison_transformers_2.dta`](../tests/data/comparison_transformers_2.dta).

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SpectrumPrecursorComparator()` | `precursor_comparator_construction_copy_and_assignment` | `name() == "SpectrumPrecursorComparator"`, one parameter, `window == ParamValue::Integer(2)` with the source's description | 4 |
| `~SpectrumPrecursorComparator()` | same | construction and drop; upstream only `delete e_ptr` | 4 |
| `SpectrumPrecursorComparator(const SpectrumPrecursorComparator&)` | same | `copy.name() == functor.name()` and equal parameters, the upstream pair | 3 |
| `SpectrumPrecursorComparator& operator=(const ...&)` | same | `clone_from` restores both | 3 |
| `double operator()(const PeakSpectrum&, const PeakSpectrum&) const` | `precursor_comparator_scores_the_parent_mass_distance`, `precursor_comparator_missing_precursor_convention_and_guards` | `1.7685` to 1e-9, `2.0` exactly, the missing-precursor and negative-window cases | 3, plus 4 |
| `double operator()(const PeakSpectrum&) const` | `precursor_comparator_scores_the_parent_mass_distance` | `self_score(a) == 2.0` exactly | 3, plus 4 |

The transcribed literal `1.7685` is tier 3. It is also derived: `DTAFile` turns
the singly-protonated masses 739.771 and 739.308 at charge 2 into
`370.3891382333855` and `370.1576382333855`, whose distance is `0.2315`, and
`2 - 0.2315 = 1.7685`. The Rust test asserts the score against that difference
computed from the loaded precursors with a tolerance of exactly zero, as well as
against the literal. Self-similarity `2.0` is derived: a zero distance leaves the
whole window, asserted with `assert_eq`.
