# BinnedSharedPeakCount support

Port of `src/openms/include/OpenMS/COMPARISON/BinnedSharedPeakCount.h` and
`src/openms/source/COMPARISON/BinnedSharedPeakCount.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `6bd5773e7fa00ea42fd8720bc85639552c3d01075b9ed927861d1d774eb09211`,
`.cpp` sha256 `ce46fcee764a2b7bd366b9e896f7f1a9f4ce3a51348362d8e614e9967471fee3`).

Rust: [`comparison::BinnedSharedPeakCount`](../src/comparison.rs).
Tests: [`tests/comparison_functors.rs`](../tests/comparison_functors.rs).
Provenance: [`tests/data/comparison_functors_provenance.json`](../tests/data/comparison_functors_provenance.json).

## API mapping

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class BinnedSharedPeakCount : public BinnedSpectrumCompareFunctor` | `pub struct BinnedSharedPeakCount` implementing [`BinnedSpectrumCompareFunctor`](BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition of a `DefaultParamHandler` instead of inheritance |
| `BinnedSharedPeakCount()` | `BinnedSharedPeakCount::new() -> Result<Self>` | fallible only because `DefaultParamHandler::new` is; the fixed name cannot be rejected. Reproduces `setName("BinnedSharedPeakCount")` then `defaultsToParam_()` (`BinnedSharedPeakCount.cpp:21-22`) |
| `BinnedSharedPeakCount(const BinnedSharedPeakCount& source)` | `Clone` | the C++ copy constructor forwards to the base and copies nothing else - notably not `precursor_mass_tolerance_` |
| `~BinnedSharedPeakCount() override` | drop glue | C++ destructor is `= default` |
| `BinnedSharedPeakCount& operator=(const BinnedSharedPeakCount& source)` | assignment of a clone | the C++ body is the self-assignment guard plus the base assignment |
| `double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override` | `BinnedSpectrumCompareFunctor::score` | see the formula and differences below |
| `double operator()(const BinnedSpectrum& spec) const override` | `BinnedSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` |
| `protected: void updateMembers_() override` | not ported | the body is empty; the trait has no `updateMembers_` hook because nothing is derived from the (empty) parameter tree |
| `protected: double precursor_mass_tolerance_` | not ported | never initialised by any constructor, never written by `updateMembers_`, never copied by the copy constructor and never read anywhere in the SDK. Reading it would be undefined behaviour; it has no observable effect and no counterpart here |

`@htmlinclude OpenMS_BinnedSharedPeakCount.parameters` documents a parameter
section that is empty, because the class registers no defaults.

## Preserved source conventions

- **Formula.** `s = spec1.bins.cwiseProduct(spec2.bins)`, then
  `static_cast<double>(s.nonZeros()) / max(spec1.nonZeros(), spec2.nonZeros())`
  (`BinnedSharedPeakCount.cpp:57-64`). The Rust port counts the indices stored in
  both maps and divides by the larger stored-bin count, converting both to `f64`
  before the division exactly as the `static_cast` does.
- **"Occupied" means stored, not nonzero.** `Eigen::SparseVector::nonZeros()`
  returns the number of *stored* coefficients, and Eigen's sparse assignment does
  not prune zeros, so a bin holding an explicit `0.0f` - written by a peak of
  intensity zero, or by a spread over an empty neighbour - counts as occupied on
  both sides of the division. `BTreeMap::len()` and key-set intersection have
  precisely that meaning.
- The `Exception::IllegalArgument` on incompatible binning
  (`BinnedSharedPeakCount.cpp:52-55`) is reproduced as `Error::InvalidValue`.
  `BinnedSpectrum::isCompatible` compares unit, size and offset, so the exception
  message "different bin size or spread" is inaccurate on both counts - spread is
  not compared, offset is - and the header's `@throw` text ("different bin size,
  offset or unit") is the correct one. The port follows the header.

## Native differences

- **Two spectra that store no bins.** The denominator is `0`, and the source
  computes `0.0 / 0` in `double`, returning NaN. The port returns `Ok(0.0)`,
  the same "defined score of 0" its sibling `BinnedSpectralContrastAngle` already
  chooses for its own degenerate case. Rule 5 of this wave forbids an unguarded
  division, and a NaN similarity silently poisons every consumer.
- Incompatible binning is `Error::InvalidValue`, not a C++ exception type.
- The comparison is refused when the two spectra together store more than
  [`MAX_COMPARED_BINS`](../src/comparison.rs) bins. The source has no ceiling.
- This score involves no floating-point accumulation at all: the numerator and
  denominator are counts, and the single division is in `f64`. It is therefore
  the one functor of the three whose result is bit-identical to the source for
  every input on which the source does not divide by zero.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| incompatible unit, size or offset | `Err(Error::InvalidValue)`, both argument orders |
| differing spread | accepted; spread is not part of compatibility, in the source or here |
| both spectra empty | `Ok(0.0)`; source returns NaN |
| one spectrum empty | `Ok(0.0)`; denominator is the other spectrum's bin count |
| stored zero bins | counted as occupied, as `nonZeros()` does |
| combined stored bins above `MAX_COMPARED_BINS` | `Err(Error::InvalidValue)` before any traversal; nothing is allocated or mutated |
| `f32` overflow | unreachable: no arithmetic on bin values |

OpenMP: no `#pragma omp` in this header or its `.cpp`.

### Class-test sections

`src/tests/class_tests/openms/source/BinnedSharedPeakCount_test.cpp`
(sha256 `4d6f983ed25bded6544315dc8758997c0ad3cee6a1f3c027a88b84337c587547`),
six sections. The fixture is `PILISSequenceDB_DFPIANGER_1.dta`, retained here as
[`tests/data/comparison_dfpianger.dta`](../tests/data/comparison_dfpianger.dta)
byte-identical (sha256
`0175e395630673e3f0a7e44ba9b450aae2af316c18d18587509fddb8e64ae4bb`).

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `BinnedSharedPeakCount()` | `binned_shared_peak_count_constructs_copies_and_assigns` | `functor.name() == "BinnedSharedPeakCount"`, parameters empty | 4 |
| `~BinnedSharedPeakCount()` | same | construction and drop; upstream only `delete ptr` | 4 |
| `BinnedSharedPeakCount(const BinnedSharedPeakCount& source)` | same | `copy.name() == functor.name()` and `copy.handler().parameters() == functor.handler().parameters()`, the upstream `TEST_EQUAL(copy.getName(), ptr->getName())` pair | 3 |
| `BinnedSharedPeakCount& operator=(const BinnedSharedPeakCount& source)` | same | assignment over a handler renamed `"scratch"` restores `"BinnedSharedPeakCount"` | 3 |
| `double operator()(const BinnedSpectrum&, const BinnedSpectrum&) const` | `binned_shared_peak_count_scores_the_upstream_fixture`, `binned_functors_reject_incompatible_binning`, `binned_shared_peak_count_counts_stored_bins_not_nonzero_values` | `score(bs1, bs2) == 0.997118` for `BinnedSpectrum(s1, 1.5, false, 2, 0)` against the same spectrum with its last peak dropped; `score(bs1, bs1) == 1.0`; a bin size of 2.0 yields `Err` | 3, plus 4 |
| `double operator()(const BinnedSpectrum&) const` | `binned_shared_peak_count_self_similarity_is_one` | `self_score(bs1) == 1.0` at offset `0.4` | 4 |

The transcribed literal `0.997118` is tier 3. It is also derived independently
and asserted exactly: the fixture bins into 347 stored bins at that layout,
dropping the last peak leaves 346, and every one of those 346 is shared, so the
score is `346.0 / 347.0 = 0.997118155...`. The Rust test asserts all three of
those counts and the exact ratio, which is stronger than the upstream
`TEST_REAL_SIMILAR`. Self-similarity `1.0` is likewise derived: the intersection
of a stored-bin set with itself is the set, so the ratio is `n/n`.
