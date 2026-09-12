# SpectrumHelper

The native `kernel::spectrum_helper` module implements the free helper
functions of `KERNEL/SpectrumHelper.h` and `KERNEL/SpectrumHelper.cpp` in Core
SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`. See the
[implementation](../src/kernel/spectrum_helper.rs),
[tests](../tests/spectrum_helper.rs) and
[hashed provenance](../tests/data/spectrum_helper_provenance.json). The source
header is template-only apart from `copySpectrumMeta`; the `.cpp` is 28 lines.
The increment adds no dependency, no feature gate and does not build or execute
C++. The source carries no `#pragma omp`; both source and port are serial.

## Public mapping

Every public member of the header is listed. `removePeaks` is mapped onto an
existing kernel operation rather than duplicated.

| Source member | Native counterpart | Status |
| --- | --- | --- |
| `template <class DataArrayT> typename DataArrayT::iterator getDataArrayByName(DataArrayT&, const std::string&)` (mutable overload) | `data_array_by_name_mut(&mut [DataArray<T>], &str) -> Option<&mut DataArray<T>>`; index form `data_array_index_by_name(&[DataArray<T>], &str) -> Option<usize>` | ported |
| `template <class DataArrayT> typename DataArrayT::const_iterator getDataArrayByName(const DataArrayT&, const std::string&)` (const overload) | `data_array_by_name(&[DataArray<T>], &str) -> Option<&DataArray<T>>`; index form `data_array_index_by_name` | ported |
| `template <typename PeakContainerT> void removePeaks(PeakContainerT& p, double pos_start, double pos_end, bool ignore_data_arrays = false)` | **mapped**: `MSSpectrum::retain_peaks` / `MSChromatogram::retain_peaks` (`src/kernel.rs`, `peak_container!` macro, `pub fn retain_peaks`) with the predicate `pos_start <= pos && pos <= pos_end`, the inclusive `PosBegin`/`PosEnd` range. Aligned arrays are kept consistent: `tests/kernel.rs::selection_and_retention_preserve_annotations_and_metadata` asserts `integer_data_arrays[0].data == vec![3]` after `retain_peaks(|peak| peak.mz > 25.0)`. `ignore_data_arrays = true` (arrays left misaligned) is not ported; see native differences | mapped with evidence |
| `template <typename PeakContainerT> void subtractMinimumIntensity(PeakContainerT& p)` | `subtract_minimum_intensity(&mut C)` and `subtract_minimum_intensity_with_limits(&mut C, SpectrumHelperLimits)` for `C: PeakContainer` | ported |
| `enum class IntensityAveragingMethod : int { MEDIAN, MEAN, SUM, MIN, MAX }` | `IntensityAveragingMethod::{Median, Mean, Sum, Min, Max}` in source order; `Default` is `Median`, the source default argument | ported |
| `template <typename PeakContainerT> void makePeakPositionUnique(PeakContainerT& p, IntensityAveragingMethod m = MEDIAN)` | `make_peak_position_unique(&mut C, IntensityAveragingMethod)` (native defaults) and `make_peak_position_unique_with(&mut C, IntensityAveragingMethod, UniquePositionOptions)`; `UniquePositionOptions::source()` reproduces the source | ported |
| `OPENMS_DLLAPI void copySpectrumMeta(const MSSpectrum& input, MSSpectrum& output, bool clear_spectrum = true)` | `copy_spectrum_meta(&MSSpectrum, &mut MSSpectrum, bool)`; the source default `true` is passed explicitly | ported |
| `PeakContainerT` template parameter (positions, intensities, three data-array lists) | `PeakContainer` trait, implemented for `MSSpectrum` (`Peak1D`, m/z) and `MSChromatogram` (`ChromatogramPeak`, RT) | native trait |
| `SpectrumHelperLimits`, `UniquePositionOptions` | native only: per-call ceilings and the lossy-behaviour opt-in | native |

The `@see makePeakPositionUnique()` on the enum is an intra-doc link. The
header's `@brief`, both `@param`s of `makePeakPositionUnique`, the `@note
Actual data is not copied` and the three `@param`s of `copySpectrumMeta` are
carried on the corresponding items. The source has no `@exception`, `@warning`
or `@deprecated` tags. The inline comment `data arrays are not updated` in
`subtractMinimumIntensity` is carried as a paragraph.

## Preserved source conventions

- `getDataArrayByName` is a linear scan by exact name; the first match wins.
  The end iterator becomes `None`.
- `subtractMinimumIntensity` takes the first minimum intensity, forms
  `rebase = -minimum` as `double`, and stores `float(intensity + rebase)`; the
  port performs the same `f32 -> f64 -> f32` arithmetic. Empty containers are
  returned unchanged. Data arrays, peak order and coordinates are untouched.
- `makePeakPositionUnique` stably sorts by position, walks the peaks once, and
  opens a new group when the position is strictly greater than the current one
  (`-0.0` and `0.0` share a group, which keeps the first position seen). Group
  intensities are widened to `double` in storage order. `Math::median` sorts
  and averages the two middle values for an even count; `Math::mean` is the
  sum divided by the count; `Math::sum` is `std::accumulate` from `0.0`;
  min/max are the extremes. The result is narrowed to the `float` peak
  intensity, as the source peak constructor does. Merged peaks come out in
  ascending position order. An empty container returns before the data-array
  policy is evaluated, exactly where the source's early return sits.
- `copySpectrumMeta` with `clear_spectrum = true` produces a spectrum holding
  only the input's metadata; with `false` the output keeps its own peaks and
  data arrays and every other field is overwritten.

## Native differences

- **Metadata survives `makePeakPositionUnique` by default.** The source ends
  with `std::swap(p_new, p)` where `p_new` is default-constructed, so the
  result loses RT, MS level, name, native ID, precursors, settings and metadata
  as well as the data arrays that the warning mentions. The port keeps every
  non-peak field. `UniquePositionOptions::reset_metadata` restores the source
  reset. This is recorded as a C++ issue candidate (the warning only announces
  the array loss).
- **Attached data arrays are refused by default.** The source logs a warning
  and drops them. Per-peak annotations cannot be merged, so the native default
  is `Error::InvalidValue`; `UniquePositionOptions::discard_data_arrays` drops
  them explicitly, and `UniquePositionOptions::source()` selects both source
  behaviours. Placeholder arrays with no entries count as attached, matching
  the source test on the array lists. No warning is logged: the crate has no
  global log stream and the refusal or the explicit option replaces it.
- **`removePeaks` is `retain_peaks`.** The source uses binary search and
  therefore requires sorted positions; `retain_peaks` is a linear predicate
  and needs no ordering. The source trims only arrays whose length equals the
  peak count and silently leaves any other array alone; `retain_peaks` rejects
  a nonempty array of the wrong length and leaves the container unchanged. The
  source `ignore_data_arrays = true` flag deliberately leaves arrays
  misaligned; that lossy option is not provided. A caller wanting it can clear
  the arrays before or after the call.
- **Finiteness is checked.** Non-finite positions or intensities are rejected
  by both mutating helpers before any write, as is a rebased or merged
  intensity that is not a finite `f32` (for example a `Sum` overflowing
  `f32`). The source performs none of these checks; a `NaN` position would
  enter `std::sort` unordered.
- **Drift time is copied.** The source `copySpectrumMeta` assigns the
  spectrum-level drift time and its unit explicitly; here
  `MSSpectrum::drift_time` and `drift_time_unit` are ordinary fields, so the
  single struct-update clone covers both. `tests/spectrum_helper.rs` asserts
  the copied value and unit. The native `peptide_identifications` field, which
  the source keeps outside `MSSpectrum`, is treated as metadata and copied.
- No `Result` is returned by `copy_spectrum_meta`; it has no failure path.

## Checked boundaries and evidence

`SpectrumHelperLimits` (`max_peaks`, `max_work`) is checked in a preflight
before any allocation or mutation. Rebasing charges `3 * n`; position merging
charges `2 * n + 8 * n * bit_length(n)`, a conservative allowance for the
stable position sort plus the per-group median sorts. Defaults are ten million
peaks and `4 * 10^9` work units, which admit the peak ceiling under the work
ceiling on 32-bit and 64-bit targets. The merged peak vector is built in a
temporary (`try_reserve_exact`) and committed only at the end, so every error
leaves the container unchanged. `data_array_*_by_name` are single linear scans
over a slice the caller already owns and take no limit. `copy_spectrum_meta`
clones metadata only, which is outside the checked-operation budgets like
ordinary `Clone`.

### Class-test sections

All nine `START_SECTION`s of `SpectrumHelper_test.cpp` are ported to
`tests/spectrum_helper.rs`; literals are transcribed (evidence tier 3, source
review). `copySpectrumMeta` has no section in the source test and receives a
native test of the `.cpp` behaviour.

| Section (source line) | Native test | Asserted values |
| --- | --- | --- |
| `MSSpectrum::FloatDataArrays getDataArrayByName` (23) | `spectrum_float_data_array_by_name` | `f1/f2/f3 -> 0/1/2`, `NOT_THERE -> None`, const overload |
| `MSSpectrum::StringDataArrays getDataArrayByName` (52) | `spectrum_string_data_array_by_name` | same indices |
| `MSSpectrum::IntegerDataArrays getDataArrayByName` (73) | `spectrum_integer_data_array_by_name` | same indices |
| `MSChromatogram::FloatDataArrays getDataArrayByName` (95) | `chromatogram_float_data_array_by_name` | same indices |
| `MSChromatogram::StringDataArrays getDataArrayByName` (117) | `chromatogram_string_data_array_by_name` | same indices |
| `MSChromatogram::IntegerDataArrays getDataArrayByName` (138) | `chromatogram_integer_data_array_by_name` | same indices |
| `removePeaks` (160) | `remove_peaks_maps_to_retain_peaks_with_inclusive_range` (through `retain_peaks`) | sizes 2/0/3/4/0, array sizes 2/0/3/4, positions 5,6 / 12,13,14 / 9..12 |
| `subtractMinimumIntensity` (211) | `subtract_minimum_intensity_source_literals` | `[-5,4] -> 0,1,9`; empty; `[5,14] -> 0,1,9` for spectrum and chromatogram |
| `makePeakPositionUnique` (257) | `make_peak_position_unique_source_literals` | MEDIAN `1,8,9,7`; MEAN `1,22/3,9,7`; SUM `1,22,9,7`; MIN `1,4,9,7`; MAX `1,10,9,7`; empty containers |

Native tests additionally cover the even-count median, signed-zero grouping,
`f64` accumulation order, metadata retention, array refusal/discard, the source
option set, non-finite and overflow rejection with bitwise unchanged checks,
limit checks before mutation, and both `copy_spectrum_meta` modes.

9 ported, 0 mapped-with-evidence, 0 mapped-without-evidence, 0 unaccounted
