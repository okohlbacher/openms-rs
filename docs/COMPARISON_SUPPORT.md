# Comparison support

This module ports selected comparison algorithms from OpenMS4-core revision `7c029e8`, except for the spectrum-similarity functor hierarchy - the two abstract bases, the three binned scorers, the alignment primitive and all seven peak-spectrum scorers - which was ported at `bc9cc12514c768385ce121d6ca4bb710fe1983c4` and carries its own support documents (`PEAK_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`, `BINNED_SPECTRUM_COMPARE_FUNCTOR_SUPPORT.md`, `BINNED_SHARED_PEAK_COUNT_SUPPORT.md`, `BINNED_SPECTRAL_CONTRAST_ANGLE_SUPPORT.md`, `BINNED_SUM_AGREEING_INTENSITIES_SUPPORT.md`, `SPECTRUM_ALIGNMENT_SUPPORT.md`, `SPECTRUM_ALIGNMENT_SCORE_SUPPORT.md`, `ZHANG_SIMILARITY_SCORE_SUPPORT.md`, `STEIN_SCOTT_IMPROVE_SCORE_SUPPORT.md`, `SPECTRUM_PRECURSOR_COMPARATOR_SUPPORT.md`, `SPECTRUM_CHEAP_DP_CORR_SUPPORT.md`, `PEAK_ALIGNMENT_SUPPORT.md`, `SPECTRAST_SIMILARITY_SCORE_SUPPORT.md`). Those scorers are now the single implementation behind every score this module offers, so the scorer and binned rows of the table below are at `bc9cc12`. It is native Rust and does not call a C++ library. Its public API uses typed configuration and fallible methods. This is algorithm-level coverage, not parity with the entire C++ comparison framework, parameter system, or all BinnedSpectrum methods.

## API mapping

| Upstream | Native Rust |
| --- | --- |
| `SpectrumAlignment` (the `DefaultParamHandler`) | `SpectrumAligner`, whose `spectrum_alignment` configures and calls `SpectrumAlignment` |
| `SpectrumAlignment::getSpectrumAlignment` | `SpectrumAlignment::align`, `Tolerance::{Absolute,Ppm}` |
| `PeakSpectrumCompareFunctor` | trait `PeakSpectrumCompareFunctor` |
| `SpectrumAlignmentScore` | `SpectrumAlignmentScore`, with `use_linear_factor` / `use_gaussian_factor` in its parameter tree |
| `KERNEL/BinnedSpectrum` | `comparison::BinnedSpectrum`, `BinConfig`, `BinUnit` |
| `BinnedSpectralContrastAngle` | `BinnedSpectralContrastAngle` (the source computes cosine, not an angle), or `binned_cosine` for the same computation without a parameter handler |
| `BinnedSharedPeakCount` | `BinnedSharedPeakCount`, or `binned_shared_peak_count` |
| `BinnedSumAgreeingIntensities` | `BinnedSumAgreeingIntensities`, or `binned_sum_agreeing_intensities` |
| `SpectrumPrecursorComparator` | `SpectrumPrecursorComparator` |
| `ZhangSimilarityScore` | `ZhangSimilarityScore` (absolute tolerance only, as supported upstream); `getFactor_` is the public `ZhangSimilarityScore::factor` |
| `SteinScottImproveScore` | `SteinScottImproveScore` |
| `SpectrumCheapDPCorr` | `SpectrumCheapDPCorr` |
| `PeakAlignment` | `PeakAlignment` |
| `SpectraSTSimilarityScore` | `SpectraSTSimilarityScore`, `SpectraStPreprocessing` |

Every one of those scorers is a `DefaultParamHandler`-carrying functor read through the `PeakSpectrumCompareFunctor` trait, exactly as the source reads `param_` inside `operator()`.

BinnedSpectrum sits in `comparison` in Rust to keep the immutable sparse representation next to its consumers. Its layout, bins, and precursor list are private and available through read-only accessors. `BinConfig::default` uses the source's recommended low-resolution width 1.0005, offset 0.4 and spread 0. For high-resolution binning use size 0.02, offset 0.0 and spread 0.

### BinnedSpectrum member review

`KERNEL/BinnedSpectrum.h` was re-reviewed member by member at Core SDK
`bc9cc12514c768385ce121d6ca4bb710fe1983c4` (header sha256
`e734cfa2a84c79d5a41e3376ad85696b9e48adb58c87f71564e5722896b5d49f`, `.cpp` sha256
`c724cc045175cdc84ff59470020461b9f74f23a96c3c6cdbb1f079348559f9e2`). The four
recommended-layout constants had no counterpart and are declared on
`BinnedSpectrum` in [`kernel/gap_closures.rs`](../src/kernel/gap_closures.rs)
(`comparison.rs` was frozen for the review; a crate-local inherent impl may
live in any module of the crate). Every public member:

| Source member | Rust counterpart | Difference |
| --- | --- | --- |
| private `MIN_MZ_ = 1.0` | literal 1.0 inside `BinConfig::bin_index` (`ln(mz)` is `ln(mz / 1.0)`) | ppm binning of m/z < 1 is a checked error; the source only asserts it in debug builds |
| `DEFAULT_BIN_WIDTH_LOWRES = 1.0005f` | `BinnedSpectrum::DEFAULT_BIN_WIDTH_LOWRES` (gap closure); also `BinConfig::default().size` | none |
| `DEFAULT_BIN_WIDTH_HIRES = 0.02f` | `BinnedSpectrum::DEFAULT_BIN_WIDTH_HIRES` | none |
| `DEFAULT_BIN_OFFSET_HIRES = 0.0f` | `BinnedSpectrum::DEFAULT_BIN_OFFSET_HIRES` | none |
| `DEFAULT_BIN_OFFSET_LOWRES = 0.4f` | `BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES`; also `BinConfig::default().offset` | none |
| `SparseVectorType = Eigen::SparseVector<float, 0, int>` | `BTreeMap<usize, f32>` | ordered map instead of a compressed sparse vector; Eigen's `nonZeros()` is `len()`, `coeffRef(i)` read is `[&i]` / `get(&i)` |
| `BinnedSpectrum()` | not ported | the source default object holds a null `bins_` pointer; `getBinIntensity` on it dereferences null. `BinnedSpectrum::new(&MSSpectrum::default(), config)` is the empty, usable equivalent |
| `BinnedSpectrum(const PeakSpectrum&, float size, bool unit_ppm, UInt spread, float offset)` | `BinnedSpectrum::new(&MSSpectrum, BinConfig { size, unit, spread, offset, .. }) -> Result` | sortedness is checked (source: debug-only precondition); bin count and spread work are bounded by `BinConfig::max_bins` / `max_updates`; f32 accumulation overflow is an error |
| copy constructor, `operator=`, `virtual ~BinnedSpectrum()` | `Clone`, assignment, drop | none |
| `operator==` (unit, size, spread, precursors, stored bins) | `PartialEq` | additionally compares `offset`, which the source omits; resource limits are ignored |
| `operator!=` | `!=` | none |
| `float getBinIntensity(double mz)` (non-const; `coeffRef` inserts a zero) | `bin_intensity(f64) -> Result<f32>` | read-only; a miss returns 0 without storing an entry; invalid m/z is an error |
| `size_t getBinIndex(float mz) const` | `config().bin_index(f64) -> Result<usize>` | narrows to f32 first, as the source signature does; negative, non-finite or out-of-range results are errors instead of a wrapped cast |
| `float getBinLowerMZ(size_t i) const` | `config().bin_lower_mz(usize) -> Result<f32>` | overflow is an error |
| `float getBinSize() const` | `config().size` | none |
| `size_t getBinSpread() const` | `config().spread` (`u32`) | none |
| `const SparseVectorType* getBins() const` | `bins() -> &BTreeMap<usize, f32>` | never null |
| `SparseVectorType* getBins()` (mutable) | not ported | bins are immutable after construction so layout and value invariants hold; the bin-wise `+` / `* 5f` algebra the class comment advertises through Eigen is not offered. Rebuild with `new` |
| `float getOffset() const` | `config().offset` | none |
| `const std::vector<Precursor>& getPrecursors() const` | `precursors() -> &[Precursor]` | none |
| `std::vector<Precursor>& getPrecursors()` (mutable) | not ported | same reason as mutable bins; construct from a spectrum with the desired precursors |
| `static bool isCompatible(const BinnedSpectrum&, const BinnedSpectrum&)` (unit, size, offset) | `a.is_compatible(&b)` | none; spread is ignored in both |
| private `binSpectrum_` | body of `new` | spread stops at bin 0 without the source's `static_cast<int>(idx - j - 1)` wraparound test |
| private `bin_spread_`, `bin_size_`, `unit_ppm_`, `offset_`, `bins_`, `precursors_` | `BinConfig` plus private `bins` and `precursors` | none |

The source class comment's `@todo` (weighted intensity spread for
high-resolution sum scores) is unimplemented upstream and here.

`BinnedSpectrum_test.cpp` (sha256
`528251907315c3d1975456c8c82dde629b53a39ebf30bb69d111e868bdba123f`) has fourteen
sections. All are ported into
[tests/binned_spectrum.rs](../tests/binned_spectrum.rs) against the retained
`comparison_dfpianger.dta` fixture (byte-identical to the source's
`PILISSequenceDB_DFPIANGER_1.dta`), tier 3 (transcribed literals `1.5`, `2`,
`1.234`, `230259`, `460519`, `690778`, `347`, `658`, `501645`, `-0.5`,
`999.5`; no C++ execution):

| Source section | Rust test |
| --- | --- |
| `~BinnedSpectrum()` | `destructor` |
| `BinnedSpectrum(const PeakSpectrum&, float, UInt, float)` | `detailed_constructor` |
| `BinnedSpectrum(const BinnedSpectrum&)` | `copy_constructor` (precursor m/z truncates to 1019) |
| `operator=(const BinnedSpectrum&)` | `assignment_operator` |
| `operator==` | `equality_operator` |
| `operator!=` | `inequality_operator` |
| `getBinSize()` | `get_bin_size` |
| `getBinSpread()` | `get_bin_spread` |
| `getBinIndex(double)` | `get_bin_index` |
| `getBinLowerMZ(size_t)` (19 assertions) | `get_bin_lower_mz` |
| `const getBins()` | `get_bins_const` |
| mutable `getBins()` | `get_bins_mutable_equivalent` (read-only lookup of bin 658 by m/z 987.0) |
| `setBinning()` (`NOT_TESTABLE`) | `set_binning_is_construction` (native spread-at-boundary check) |
| `isCompatible(a, b)` | `is_compatible` |
| (none) | extra: `recommended_layout_constants` |

Self-audit (`BinnedSpectrum.h`): 14 ported, 0 mapped-with-evidence, 0
mapped-without-evidence, 0 unaccounted. The earlier golden and boundary tests in
[tests/comparison.rs](../tests/comparison.rs) remain and are not counted here.

## Alignment and score semantics

Absolute alignment preserves the actual source band traversal, fallback costs for absent matrix cells, diagonal-before-gap tie ordering, and traceback starting at the last selected diagonal. It is **not** a greedy nearest-peak substitution. Returned pairs are ordered, one-to-one, and within the inclusive tolerance. Compact contiguous row storage replaces C++ nested maps; 1,000 deterministic cases compare this representation against a separate map-shaped source transcription.

Ppm alignment preserves `MatchedIterator` traversal, including its float32 distance/tolerance calculations, lower-position tie preference and target reuse. The result is directed: exchanging reference and target can change the matches. The source stops when a target distance does not strictly improve; duplicate targets can therefore prevent it from reaching a later closer target. This behavior is explicitly tested and retained instead of claiming unconditional nearest-neighbor optimality. Ppm tolerance windows are computed from the reference m/z.

AlignmentScore uses `sum(sqrt(I1 * I2 * factor)) / sqrt(sum(I1²) * sum(I2²))`, with the source's grouping: **one** `sqrt` of the product of the two squared-intensity sums, not the product of two `sqrt`s. A self-score need not be one; scores can exceed one. Linear weighting is `(tolerance-distance)/tolerance`; Gaussian weighting is `erfc(distance/(3*tolerance*sqrt(2)))`, evaluated by `libm::erfc` as `std::erfc`.

A **zero tolerance under a weighting flag** is an error, not a factor of one. `diff_align <= tolerance` in `SpectrumAlignment.h` is inclusive, so two peaks at the same m/z are still aligned, and both weighting expressions then divide by zero: the linear factor is `(0.0 - 0.0) / 0.0` and the Gaussian argument is `0.0 / (3.0 * 0.0 * sqrt(2))`. A statement-for-statement C++ transcription of `SpectrumAlignment.h` and `SpectrumAlignmentScore.cpp`, run on `s1 == s2 == {(100.0, 1.0f)}` at tolerance `0`, prints `factor = nan` and `score = nan` under either flag and `score = 1` under neither. Returning `1` there - which an earlier port did, by special-casing a zero distance inside a zero window - reports the *unweighted* score under a weighting configuration; the port now refuses. Unweighted, the same input still scores `1`. If ppm float rounding admits a match outside the exact weighting window, weighted scoring likewise returns an error instead of a negative radicand's NaN; the Gaussian arm stays finite there and scores below one.

Zhang includes **all** pairs at distance strictly smaller than its absolute tolerance and normalizes by `sqrt(TIC1 * TIC2)`. Stein/Scott includes all pairs at distance **at most twice** its tolerance, subtracts `z = (tolerance/10000) * (TIC1 * TIC2)` with that grouping, normalizes by `sqrt(sum(I1²) * sum(I2²))`, and applies the configured threshold after narrowing it to `float` as the source's `(float)param_.getValue("threshold")` does. Neither is replaced with a one-to-one match score. Both walk the source's sticky `j_left` cursor, which is read when a reference peak's inner loop starts and written during it, so a write takes effect on the *next* reference peak.

The precursor score is `max(0, window - abs(mz1-mz2))`, using the first precursor. As in the source, an absent precursor is treated as m/z zero, so two missing precursors score the full window. Charge does not participate in this score.

## Binning and numerical policy

Absolute bins use `floor(float32(mz)/size + offset)` with source float32 arithmetic. Ppm bins use the logarithmic source convention with minimum m/z one; the numerator preserves the float overload of `log`. Spread copies each intensity to neighboring bins and stops at the left boundary, so spreading increases the sum of stored intensities. Bin accumulation is float32, with overflow returned as an error.

Explicit zero coefficients remain stored, including cancellation to zero. Shared-bin count uses the intersection of **stored indices**, divided by the larger stored-index count, matching Eigen's sparse-storage behavior; it is not a count of strictly positive peaks. A lookup does not insert a zero coefficient in Rust. Layout compatibility compares size, unit and offset and ignores spread, as in the source. Equality includes offset (fixing its omission in the source equality operator), includes spread/data/precursors, and ignores resource budgets.

Peak-level score products reproduce the source exactly: `Peak1D::getIntensity` returns `float`, so `s1[i].getIntensity() * s2[j].getIntensity()` is a single-precision multiply whose *result* is widened when it meets a `double`, and every reduction is a sequential float64 accumulation in source order. Computing that product in float64 would be exact and would therefore **not** match the source; it is worth about 3.5e-9 relative on the `DFPIANGER` fixture. The three **binned** scores reproduce the source's float32 reduction for the same reason, because `Eigen::SparseVector<float>::dot` and `::sum` return `float` and the source's `const double` results are widenings of that; see the three `BINNED_*_SUPPORT.md` documents and the reduction-association caveat below. The current source fixture gives alignment self-score `1.484501010820065`, independently checked from the formula; its historical upstream literal `1.48268` passes that test's 0.01 absolute tolerance. Tests retain this loose historical comparison and add the tight independent expectation. The truncated-reference alignment score is `3.8247227487373805`. Computing the same two numbers with a float64 intensity product gives `1.4845010143546342`, which the class test's own 0.01 band cannot tell apart - the class test alone does not pin this down, which is why the tight literal is asserted. Binned golden comparisons likewise permit small differences from float32 reduction order.

All alignment/scoring inputs require sorted, finite, nonnegative m/z and consistent auxiliary array lengths. Alignment and binning permit signed finite intensities. The binned contrast angle supports signed bins and does not clamp, as the source does not. Square-root peak scores reject negative intensities; the binned agreeing-intensity score accepts negative bins and truncates them away, as the source does. Empty or zero-norm scores return zero instead of NaN; a stored all-zero binned spectrum can still have shared-bin self-score one. Comparisons never mutate inputs.

**One implementation per peak-spectrum score.** An earlier wave shipped `SpectrumAlignmentScore`, `ZhangSimilarityScore`, `SteinScottImproveScore` and `SpectrumPrecursorComparator` as typed `Copy` configuration structs with no `DefaultParamHandler`, and a later wave added faithful `DefaultParamHandler` functors beside them under `-er` names. The two disagreed: the struct divided by `sqrt(sum1) * sqrt(sum2)` where `SpectrumAlignmentScore.cpp:99` and `SteinScottImproveScore.cpp:113` write `sqrt(sum1 * sum2)`, formed the intensity product in `f64` where the source multiplies two `float`s, re-derived the many-to-many cursor per reference peak instead of carrying the source's sticky `j_left`, and grouped Stein/Scott's subtracted term differently from `constant * (sum3 * sum4)`. The structs are gone. Each of these headers now has exactly one implementation in this crate, it is the faithful functor, and it carries the header's own name.

**One implementation per binned score.** `binned_cosine`, `binned_shared_peak_count` and `binned_sum_agreeing_intensities` used to be separate float64 ports that disagreed with the `BinnedSpectralContrastAngle` / `BinnedSharedPeakCount` / `BinnedSumAgreeingIntensities` functors in reduction precision, denominator grouping, clamping and degenerate-case policy. Each pair is now a single implementation with two entry points: the functor is the port of the header and is authoritative, and the function is the same code path for callers that want neither a parameter handler nor a `&dyn BinnedSpectrumCompareFunctor`. `tests/comparison_functors.rs::the_parameter_free_functions_are_the_functor_implementations` asserts the agreement across fixture, degenerate, negative-bin and incompatible inputs, and `the_surviving_contrast_angle_policy_is_the_source_one` asserts that the retained float32 result genuinely differs from the float64 one it replaced.

**Reduction association.** The port sums float32 coefficients sequentially in ascending bin index. That is exactly what Eigen's sparse `dot` does (a scalar merge of two inner iterators), but *not* what `SparseVector::sum()` and `s.coeffs().cwiseMax(0).sum()` do: both map the stored-value array to a dense vector and use the vectorised dense reduction, whose association depends on packet width and Eigen version. For those two reductions - the totals and the numerator of `BinnedSumAgreeingIntensities` - what the port reproduces is the float32 *precision*, and agreement with a vectorised C++ build is to float32 reduction rounding rather than bit for bit. A C++ build that contracts `res += a * b` into an FMA (the default under `-ffp-contract=fast`) can likewise differ in the last bit of the dot product. `tests/comparison_functors.rs::binned_reductions_are_sequential_f32_in_ascending_bin_order` pins the association the port uses, which is also the one a parallel reduction would have to reproduce under `src/concept/parallel.rs`.

The Zhang Gaussian width is instance-specific. This deliberately corrects the C++ function-local static denominator that caches the tolerance of its first call.

## Resource limits

`SpectrumAlignment::max_cells` defaults to 5,000,000 and counts DP cells plus row/column initialization (or the ppm input size). Worst-case absolute alignment remains O(n*m) time and stored cells; its band often reduces both. Ppm traversal is O(n+m). `BinConfig` defaults to 1,000,000 stored bins and 10,000,000 conservative spread updates. Index and arithmetic overflow are checked. Sparse storage avoids allocating up to the highest bin index. Zhang and Stein/Scott default to 5,000,000 candidate pairs; rejected tolerance-boundary candidates count toward this work limit as well.

## Validation and provenance

`cargo test --offline --no-default-features --test comparison` covers source alignment fixtures, float ppm ties/duplicates, the alignment/score distinction, binned counts and values, precursor scores, Zhang/SteinScott golden values, Gaussian factors, the zero-tolerance weighting NaN, empty/zero/signed inputs, compatibility, resource limits, and the independent DP oracle. `--test comparison_scorers` covers every `START_SECTION` of the eight upstream class tests behind the alignment primitive and the seven scorers. `cargo clippy --offline --no-default-features --lib --test comparison -- -D warnings` passes. No C++ build was run.

The four `comparison_*.dta` fixtures are unmodified copies of `SpectrumAlignment_in1.dta`, `SpectrumAlignment_in2.dta`, `PILISSequenceDB_DFPIANGER_1.dta`, and `Transformers_tests_2.dta` from upstream `src/tests/class_tests/openms/data/`. They retain the upstream BSD-3-Clause project license. Tests also reuse the previously copied `Transformers_tests.dta`.

| Local fixture | SHA-256 |
| --- | --- |
| `comparison_alignment_1.dta` | `10986501b97f6c6de9a37808c97134c9572a5a179465e7034795cb724395933e` |
| `comparison_alignment_2.dta` | `9aea8ed6847a47918ecc906b156cf50f75605a02d4284cd1de454e1845090a3f` |
| `comparison_dfpianger.dta` | `0175e395630673e3f0a7e44ba9b450aae2af316c18d18587509fddb8e64ae4bb` |
| `comparison_transformers_2.dta` | `a21f79a751c8ac69a2c3e35a26da2f77bfd9647a7baefc85a4ee5516257adae5` |

| Upstream source path | SHA-256 |
| --- | --- |
| `src/openms/include/OpenMS/COMPARISON/SpectrumAlignment.h` | `e38817397c8f34619b0231887d0b493928f8889163da88ab6eee275e8313e021` |
| `src/openms/include/OpenMS/DATASTRUCTURES/MatchedIterator.h` | `79941c1975151fcb2b98228295fcdbad0c8aca966266fc20b053ee27aeeab0d3` |
| `src/openms/source/COMPARISON/SpectrumAlignmentScore.cpp` | `2e08a1e8650f3b071068ab0094682ad5b01f2ac0dddcad8cc8c98d97cf8ea111` |
| `src/openms/source/KERNEL/BinnedSpectrum.cpp` | `c724cc045175cdc84ff59470020461b9f74f23a96c3c6cdbb1f079348559f9e2` |
| `src/openms/include/OpenMS/KERNEL/BinnedSpectrum.h` | `e734cfa2a84c79d5a41e3376ad85696b9e48adb58c87f71564e5722896b5d49f` |
| `src/openms/source/COMPARISON/BinnedSpectralContrastAngle.cpp` | `ea854a8917eedb06c6fa44bc0b49bea6673cac53da3f6f260409f1af34d66ac6` |
| `src/openms/source/COMPARISON/BinnedSharedPeakCount.cpp` | `ce46fcee764a2b7bd366b9e896f7f1a9f4ce3a51348362d8e614e9967471fee3` |
| `src/openms/source/COMPARISON/BinnedSumAgreeingIntensities.cpp` | `a0f3e2553b3858c2d8dc8cfa5340c2221b208d989c1b7978988ba349bb67fff0` |
| `src/openms/source/COMPARISON/SpectrumPrecursorComparator.cpp` | `eec1ce3b572a00631f3198c1c6717928c136cc93a2c4d1aff3621508d32b2527` |
| `src/openms/source/COMPARISON/ZhangSimilarityScore.cpp` | `66aa29df1e5a4ce270a8b2e69aa2849b190fe6071c71f2c4a05b70183929c7d3` |
| `src/openms/source/COMPARISON/SteinScottImproveScore.cpp` | `a6f2b273987ef23d86a12f70cfb2c45cf3ad7c616b57a1d623af07684bfeb84d` |
| `src/openms/source/COMPARISON/SpectrumAlignment.cpp` | `34c39c22df3d5cb7d5d2781a39959023d80084277c40db9431d6ca8a4f3af008` |
| `src/openms/source/COMPARISON/SpectrumCheapDPCorr.cpp` | `5bad4274f5b7bc69cbc9e8696535f4e718529d987367b9695760a839ea2bff4f` |
| `src/openms/source/COMPARISON/PeakAlignment.cpp` | `a621dee886979eba7771ddbe52cc8c450d108e872971baef395542e04979a053` |
| `src/openms/source/COMPARISON/SpectraSTSimilarityScore.cpp` | `f017f2adcfeab5edf829641066ad636cd6c378e397d55531f1588126aa7f93a5` |
| `src/openms/include/OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h` | `3e5d1af32c7704e4368398ac5742603ac876405494bb13899216a9f6c1c77fcc` |
