# SpectraSTSimilarityScore support

Port of `src/openms/include/OpenMS/COMPARISON/SpectraSTSimilarityScore.h` and
`src/openms/source/COMPARISON/SpectraSTSimilarityScore.cpp` at Core SDK revision
`bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(header sha256 `131d8833334424c80faaf1d2f33569f259f9bfe7cd29773161d0826cc533b727`,
`.cpp` sha256 `f017f2adcfeab5edf829641066ad636cd6c378e397d55531f1588126aa7f93a5`).

Rust: [`comparison::SpectraSTSimilarityScore`](../src/comparison.rs).
Tests: [`tests/comparison_scorers.rs`](../tests/comparison_scorers.rs).
Provenance: [`tests/data/comparison_scorers_provenance.json`](../tests/data/comparison_scorers_provenance.json).

## API mapping

Every public member of the header appears here.

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| `class SpectraSTSimilarityScore : public PeakSpectrumCompareFunctor` | `pub struct SpectraSTSimilarityScore` implementing [`PeakSpectrumCompareFunctor`](PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md) | composition of a `DefaultParamHandler` instead of inheritance |
| `SpectraSTSimilarityScore()` | `SpectraSTSimilarityScore::new() -> Result<Self>` | reproduces `setName("SpectraSTSimilarityScore")` **without** the `defaultsToParam_()` every sibling calls (`SpectraSTSimilarityScore.cpp:18-22`); harmless only because the class registers no defaults |
| `SpectraSTSimilarityScore(const SpectraSTSimilarityScore& source)` | `Clone` | the C++ copy constructor is `= default` |
| `~SpectraSTSimilarityScore() override` | drop glue | C++ destructor is `= default` |
| `SpectraSTSimilarityScore& operator=(const SpectraSTSimilarityScore& source)` | `Clone::clone_from`, or assignment of a clone | the C++ body is the self-assignment guard plus the base assignment |
| `double operator()(const PeakSpectrum& spec1, const PeakSpectrum& spec2) const override` | `PeakSpectrumCompareFunctor::score` | exactly `dot(transform(a), transform(b))`; the source spells the binning and normalisation out a second time with the same arguments |
| `double operator()(const BinnedSpectrum& bin1, const BinnedSpectrum& bin2) const` | `dot(&self, &BinnedSpectrum, &BinnedSpectrum) -> Result<f64>` | a distinct name, because Rust has no overloading; adds the compatibility and size checks the source omits |
| `double operator()(const PeakSpectrum& spec) const override` | `PeakSpectrumCompareFunctor::self_score`, the trait default | the C++ override is `return operator()(spec, spec)` (`SpectraSTSimilarityScore.cpp:37-40`) |
| `bool preprocess(PeakSpectrum& spec, float remove_peak_intensity_threshold = 2.01, UInt cut_peaks_below = 1000, Size min_peak_number = 5, Size max_peak_number = 150)` | `preprocess(&self, &mut MSSpectrum, SpectraStPreprocessing) -> Result<bool>` | the four defaulted arguments become one `Copy` options struct with the same defaults, because Rust has no default arguments |
| `BinnedSpectrum transform(const PeakSpectrum& spec)` | `transform(&self, &MSSpectrum) -> Result<BinnedSpectrum>` | returns `Result`; a zero norm is refused instead of filling the vector with NaN |
| `double dot_bias(const BinnedSpectrum& bin1, const BinnedSpectrum& bin2, double dot_product = -1) const` | `dot_bias(&self, &BinnedSpectrum, &BinnedSpectrum, Option<f64>) -> Result<f64>` | the `-1` sentinel becomes `Option`; any non-positive `Some` still recomputes, exactly as the source's `(dot_product > 0) ? ... : ...` does |
| `double delta_D(double top_hit, double runner_up)` | `delta_d(&self, f64, f64) -> Result<f64>` | `Exception::DivisionByZero` becomes `Error::InvalidValue`; the source's method is non-const and touches no state, so this takes `&self` |
| `double compute_F(double dot_product, double delta_D, double dot_bias)` | `compute_f(&self, f64, f64, f64) -> Result<f64>` | returns `Result` so that a non-finite argument is refused rather than silently selecting `b = 0` |
| the fixed binning `BinnedSpectrum(spec, 1, false, 1, DEFAULT_BIN_OFFSET_LOWRES)` | `SpectraSTSimilarityScore::bin_config() -> BinConfig` | the source repeats the five arguments at `.cpp:45`, `:46` and `:96`; here they are named once |
| `protected:` (empty) | nothing | the header's `protected` section declares no members |

The class registers **no parameters**, so there is no `@htmlinclude` block and
the handler carries only the name.

## The scaling exponents

SpectraST scales intensity by `0.5` and m/z by `0`, and this implementation does
exactly that:

* `preprocess` replaces every surviving intensity with `sqrt(intensity)`,
  computed in `float` because `Peak1D::IntensityType` is `float` and the C++
  overload resolves to `std::sqrt(float)` (`SpectraSTSimilarityScore.cpp:81`).
* **No mass weighting is applied anywhere.** The m/z of a surviving peak is
  copied unchanged, so the m/z exponent is `0`. The `m/z^0.5` variant that
  appears in descriptions of the published score is not implemented upstream, and
  introducing it here would change every score this class has ever produced.

The only other transform is the normalisation in `transform`, which divides the
binned vector by its own Euclidean norm so that a spectrum scores exactly `1`
against itself up to the `f32` reduction's rounding.

## Preserved source conventions

- **Binning is fixed**: width `1` Th, absolute units, spread `1`, offset
  `BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES = 0.4f`, under the source's own
  `// TODO: resolution seems rather low`. A peak therefore lands in
  `floor(mz + 0.4)` and also contributes to the neighbouring bin on each side.
- **Eigen's reductions in `f32`.** `dot` is Eigen's sparse dot: the products of
  the coefficients stored at indices present in both vectors, accumulated in
  `f32` in ascending index order and widened only on return - the same helper the
  wave-A binned scorers use. `norm()` is `sqrt(cwiseAbs2().sum())` with the
  squares and their accumulation in `f32` too, and `operator/=` divides each
  stored coefficient in place, so stored zeros stay stored.
- **`dot_bias`'s sentinel.** `(dot_product > 0) ? dot_product : (*this)(bin1, bin2)`
  treats the documented `-1` and the legacy `0` alike as "recompute", and the
  following `if (denominator <= 0.0) return 0.0` is the source's own guard
  against a division by zero. Both are reproduced.
- **`compute_F`'s bands.** `b` is `0.12` below `0.1` *and* on `(0.35, 0.4]`,
  `0.18` on `(0.4, 0.45]`, `0.24` above `0.45`, and zero on `[0.1, 0.35]`. The
  low-bias and high-bias penalties sharing a value is what
  `SpectraSTSimilarityScore.cpp:132` writes. The final expression is
  `0.6 * dot_product + 0.4 * delta_D - b`.
- **`preprocess` keeps an m/z prefix, not the strongest peaks.** The header says
  it "cuts peaks exceeding the max_peak_number most intense peaks"; the code
  sorts by **position** and stops after examining `max_peak_number` peaks
  (`SpectraSTSimilarityScore.cpp:72-87`), so `max_peak_number` caps how many
  peaks are looked at and the result can be shorter. Reproduced, and documented
  at the item.
- **`preprocess` replaces the whole spectrum.** `spec = tmp` assigns a
  default-constructed `PeakSpectrum` over the argument, discarding precursors,
  retention time, native id and every data array. This is the one deliberately
  preserved lossy behaviour in the package: a SpectraST workflow's numbers depend
  on which peaks survive, and quietly keeping more state would be a different
  function. The rustdoc says so in a dedicated section, and a test pins it.
- **The relative intensity floor** is `(1.0 / cut_peaks_below) * max_intensity`,
  with the base peak found the way `std::max_element` with `Peak1D::IntensityLess`
  finds it - the first maximal peak.

## Native differences

- `dot` refuses incompatible binning and charges the pair against
  `MAX_COMPARED_BINS` before traversing anything. The source hands mismatched
  vectors straight to Eigen, which asserts in a debug build and reads past the
  shorter one otherwise.
- `transform` refuses a zero or non-finite norm - an empty spectrum, or one whose
  intensities are all zero. The source divides by that zero and fills the vector
  with NaN, which then poisons every score computed from it.
- `preprocess` refuses a `cut_peaks_below` of zero, where the source's
  `1.0 / cut_peaks_below` is `inf` and every peak is discarded, and refuses a
  negative m/z or intensity, where `sqrt` would yield NaN. The spectrum is left
  untouched when any check fires.
- `delta_d` refuses a non-finite argument as well as the source's zero `top_hit`.
- `compute_f` refuses non-finite arguments. The source performs no check and lets
  a NaN `dot_bias` select `b = 0`, because every comparison against NaN is false.
- `dot_bias` refuses a non-finite supplied `dot_product`, which slips past the
  source's `<= 0.0` guard and produces NaN.
- The `f32` accumulations are checked for overflow rather than carrying an
  infinity into the score.

### Reduction fidelity

`dot` is Eigen's scalar sparse-dot merge and is reproduced exactly. `norm()` and
the `dot_bias` numerator reduce in ascending bin order, which reproduces the
`f32` *precision* of Eigen's vectorised dense redux over the stored value array
rather than its bit pattern on every build - the same bounded claim
`BINNED_SUM_AGREEING_INTENSITIES_SUPPORT.md` makes, for the same reason. A C++
build contracting `res += a * b` into an FMA can also differ in the last bit.

## Checked boundaries and evidence

| Boundary | Behaviour |
| --- | --- |
| incompatible binning in `dot` or `dot_bias` | `Err(Error::InvalidValue)`; source hands it to Eigen |
| combined stored bins above `MAX_COMPARED_BINS` | `Err(Error::InvalidValue)` before traversal |
| `transform` of an empty or all-zero spectrum | `Err(Error::InvalidValue)`; source yields NaN |
| unsorted peaks | `Err(Error::UnsortedData)` from the binning |
| disjoint bins | `Ok(0.0)` exactly; the merge finds no shared index |
| `dot_bias` with a non-positive `dot_product` | recomputes it, then returns `0.0` if that is still not positive, as upstream |
| `dot_bias` with a NaN `dot_product` | `Err(Error::InvalidValue)`; source returns NaN |
| `delta_d(0, x)` | `Err(Error::InvalidValue)`; source throws `DivisionByZero` |
| `compute_f` with a NaN term | `Err(Error::InvalidValue)`; source selects `b = 0` |
| `preprocess` with `cut_peaks_below == 0` | `Err(Error::InvalidValue)`, spectrum untouched; source discards every peak |
| `preprocess` metadata | discarded, as upstream; documented at the item and tested |
| `f32` overflow in any reduction | `Err(Error::InvalidValue)` |

OpenMP: no `#pragma omp` in this header or its `.cpp`; the port is serial, as the
source is.

### Class-test sections

`src/tests/class_tests/openms/source/SpectraSTSimilarityScore_test.cpp`
(sha256 `fff8b7cde10809ee390d4d8673cc57bd8e07056b791fb8deb28876d91965c8f0`),
**twelve sections, all ported**, one of them `NOT_TESTABLE` upstream. The fixture
is `SpectraSTSimilarityScore_1.msp`, retained byte-identical as
[`tests/data/comparison_spectrast_1.msp`](../tests/data/comparison_spectrast_1.msp)
(sha256 `3136cf6ab0eda25d2ade5d746a76b33a326662caf19dec69ed8c7efea6024ebb`).
`MSPFile` is not ported, so the test reads only the `Num peaks:` blocks and the
tab-separated peak lines that follow them; that is a fixture reader, not an MSP
format port, and it is stated as such in the test.

| Section | Rust test | Asserted value | Tier |
| --- | --- | --- | --- |
| `SpectraSTSimilarityScore()` | `spectrast_construction_copy_and_assignment` | `name() == "SpectraSTSimilarityScore"`, empty parameters and empty defaults | 4 |
| `~SpectraSTSimilarityScore()` | same | construction and drop | 4 |
| `SpectraSTSimilarityScore(const SpectraSTSimilarityScore&)` | same | equal name and parameters, the upstream pair | 3 |
| `SpectraSTSimilarityScore& operator=(const ...&)` | same | equal name and parameters | 3 |
| `double operator()(const PeakSpectrum&) const` | `spectrast_peak_spectrum_dot_products` | `self_score(s1) == 1` to 1e-6 | 3, plus 4 |
| `double operator()(const PeakSpectrum&, const PeakSpectrum&) const` | same | `score(s1, s2) == 1` to 1e-6 for the two identical MSP spectra, `score(s1, s3) == 0` exactly | 3, plus 4 |
| `(double operator()(const BinnedSpectrum&, const BinnedSpectrum&) const)` | same | the same two values through `dot(transform, transform)`, their equality with the peak-spectrum overload, and the incompatible-binning refusal | 3, plus 4 |
| `bool preprocess(...)` | `spectrast_preprocess_filters_and_square_roots_the_intensities` | `6` survivors at threshold 2, `false` at `min_peak_number` 12, `8` at `max_peak_number` 8, the `sqrt` intensity, the m/z-prefix behaviour, the discarded metadata and the zero-`cut_peaks_below` refusal | 3, plus 4 |
| `double delta_D(double, double)` | `spectrast_delta_d_and_compute_f` | `0.2`, `0.96`, and the zero and NaN refusals | 3, plus 4 |
| `(double compute_F(double, double, double))` | same | upstream is `NOT_TESTABLE`; seven points pin all four bias bands and both band edges | 4 |
| `double dot_bias(const BinnedSpectrum&, const BinnedSpectrum&, double)` | `spectrast_dot_bias_measures_domination_by_few_bins` | `98.585`, symmetry, the orthogonal `0.0`, both sentinels, a halved denominator and the NaN refusal | 3, plus 4 |
| `BinnedSpectrum transform(const PeakSpectrum&)` | `spectrast_transform_normalises_the_binned_vector` | `0.1205`, `0.3614`, `0.602`, `0.602`, and every coefficient exactly against `raw / sqrt(69)` | 3, plus 4 |

The transcribed literals are tier 3 and every one of them is also derived:

* `transform`: `floor(mz + 0.4)` puts the four peaks of the test spectrum in bins
  `0, 1, 2, 3`, and spread `1` adds each intensity to both neighbours with bin 0
  unable to spread downwards, giving stored bins `1, 3, 5, 5, 3` and a norm of
  `sqrt(69)`. The Rust test asserts the bin indices, the raw values and the exact
  quotients, which is stronger than the upstream `TEST_REAL_SIMILAR`.
* `dot_bias`: the two raw binned vectors are `(1, 1, 3, 5, 5, 3)` and
  `(0, 4, 9, 15, 11, 6, 0)`; their shared-index products are `0, 4, 27, 75, 55,
  18` and the sum of their squares is `9719`, so the value is the `f32`
  `sqrt(9719)` - asserted exactly, not to a tolerance.
* `score(s1, s3) == 0`: with bin width 1 and spread 1 the two peak sets occupy
  disjoint bins, so the intersection the sparse dot walks is empty and the result
  is exactly zero.
* `preprocess`: the ten MSP intensities are `2, 2, 5, 2, 3, 4, 5, 5, 2, 3` and
  the relative floor at `cut_peaks_below = 10000` is `5e-4`, so exactly the six
  peaks above intensity 2 survive; at `max_peak_number = 8` all eight examined
  peaks pass.
* `delta_D`: `(5 - 4) / 5` and `(25 - 1) / 25`.
