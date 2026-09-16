# Signal-to-noise estimation: the median estimator and its base

The reference is OpenMS4-core commit `bc9cc12514c768385ce121d6ca4bb710fe1983c4`
(`.reference/openms4-core-bc9cc12/`). The reference platform for numerics is
the Linux x86-64 **Release** build `openms4-release-bc9cc12-c19e494-174b576`
(conda-forge GCC 14.4.0). No C++ is built or called by the crate.

| Source | Rust |
|---|---|
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h` (a header template) and `source/.../SignalToNoiseEstimatorMedian.cpp` | `src/processing/peak_picking/noise.rs` |
| `PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h` | `src/processing/noise_estimation.rs` (and the `estimate_` use in `noise.rs`) |
| `source/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.cpp` (`estimateNoiseFromRandomScans`) | `src/processing/noise_estimation.rs` |

`SignalToNoiseEstimatorMeanIterative.h`, the other subclass in the pinned core,
is ported in `processing::mean_noise` (`MEAN_NOISE_SUPPORT.md`); this document
covers only the trait it implements. The picker that consumes the median
estimator is `PEAK_PICKING_SUPPORT.md`.

```rust
use openms::processing::peak_picking::{PickingCompatibility, SignalToNoiseEstimatorMedian};

fn ratios(positions: &[f64], intensities: &[f64]) -> openms::Result<Vec<f64>> {
    let estimator = SignalToNoiseEstimatorMedian {
        window_length: 40.0,
        min_required_elements: 10,
        ..Default::default()
    };
    // The TOPP tools use the source profile; `estimate` is the native one.
    let estimates = estimator.estimate_with_compatibility(
        positions,
        intensities,
        &PickingCompatibility::source(),
    )?;
    Ok(estimates.signal_to_noise)
}
```

## Two profiles

`PickingCompatibility::default()` is the **native safety profile**: it refuses
inputs and parameter values on which the source computes non-finite or
degenerate results and bins an out-of-range intensity by the documented intent.
`PickingCompatibility::source()`, which every TOPP tool path uses
(`src/cli/tools/peak_picker_hi_res.rs:415`), computes what the Release build
computes wherever the source is defined and refuses exactly where it is
undefined. The estimator's own switches are `NoiseCompatibility`
(`source_value_domain`, `nan_for_empty_input`, `bin_index`), a field of
`PickingCompatibility`, in the pattern of the FFAP module's
`AbundanceOverride`.

## API mapping

### `SignalToNoiseEstimatorMedian.h`: the 13 public declarations

| # | Source (line) | Rust | Notes |
|---|---|---|---|
| 1 | `enum IntensityThresholdCalculation {MANUAL = -1, AUTOMAXBYSTDEV = 0, AUTOMAXBYPERCENT = 1}` (:68) | `NoiseHistogramRange::{Manual, StandardDeviation, Percentile}` | The variant carries the value its mode reads. `AUTOMAXBYPERCENT` is ported on its defined domain (below). |
| 2 | `using SignalToNoiseEstimator<Container>::stn_estimates_` (:70) | `NoiseEstimates::signal_to_noise`, `SignalToNoiseEstimator::signal_to_noise` | The port is stateless; see *State* under [Native differences](#native-differences). |
| 3 | `using ...::defaults_` (:71) | `SignalToNoiseEstimatorMedian::defaults()`; `DefaultParamHandler::defaults` of `from_param_with_handler` | Equal trees. |
| 4 | `using ...::param_` (:72) | `DefaultParamHandler::parameters` of `from_param_with_handler`; `to_param` | The handler holds the caller's tree merged with the defaults, as `param_`; `to_param` derives the default layout from the typed fields. Values agree; tags, descriptions and restrictions the caller changed are kept only by the handler. As in the source, editing the tree does not change the estimator. |
| 5 | `typedef ...::PeakIterator PeakIterator` (:74) | the entry points `estimate` (slices), `estimate_spectrum`, `estimate_chromatogram`, `estimate_peaks<P: NoisePoint>` | `Container::const_iterator`; Rust takes slices. |
| 6 | `typedef ...::PeakType PeakType` (:75) | trait `NoisePoint` (`position`, `intensity`), implemented for `Peak1D` and `ChromatogramPeak` | Reduced to the two members `computeSTN_` reads. |
| 7 | `typedef ...::GaussianEstimate GaussianEstimate` (:77) | `noise_estimation::GaussianEstimate` | Public here; protected struct in the base. |
| 8 | `SignalToNoiseEstimatorMedian()` (:80) | `SignalToNoiseEstimatorMedian::default()` | Also sets the name, `SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME`. |
| 9 | copy constructor (:123) | `Clone` | The source copy clears `stn_estimates_` (`updateMembers_`); a clone has none. |
| 10 | `operator=` (:133) | `Clone::clone_from` | As 9; the percentages are not assigned in the source (see State). |
| 11 | `~SignalToNoiseEstimatorMedian()` (:145) | `Drop` (implicit) | |
| 12 | `getSparseWindowPercent()` (:149) | `NoiseEstimates::sparse_window_percent` | |
| 13 | `getHistogramRightmostPercent()` (:155) | `NoiseEstimates::histogram_rightmost_percent` | |

The count is 13: the enum, three `using` re-exports, three typedefs, four
special members and two getters (the earlier table mapped 8 and missed rows 2 to
6).

### `SignalToNoiseEstimatorMedian.h`: protected members

| Source | Rust | Notes |
|---|---|---|
| `computeSTN_(const Container&)` (:168) | private `estimate_points`, `windows`, `percentile_upper_end`; public entry points above | |
| `updateMembers_()` (:398) | `from_param`, `from_param_with_warnings`, `from_param_with_handler`, `set_parameters` (private `from_complete_param`) | `auto_mode` 0 and 1 select their modes; every other value is manual, and only `-1` passes the restriction. |
| `max_intensity_`, `auto_max_stdev_Factor_`, `auto_max_percentile_`, `auto_mode_` | `histogram_range` plus `range_parameters: NoiseRangeParameters` | The record keeps the values the selected mode ignores. The value after estimation is `NoiseEstimates::max_intensity`. |
| `win_len_`, `bin_count_`, `min_required_elements_`, `noise_for_empty_window_`, `write_log_messages_` | `window_length`, `bin_count`, `min_required_elements`, `noise_for_empty_window`, `write_log_messages` | `write_log_messages` gates two of the three warnings, as in the source. |
| `sparse_window_percent_`, `histogram_oob_percent_` (:434-436) | `NoiseEstimates` | |

### Inherited members the estimator uses

| Source | Rust |
|---|---|
| `DefaultParamHandler::setName`, `getDefaults`, `setParameters`, `getParameters`, `defaultsToParam_` | `SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME`, `defaults`, `from_param*`/`set_parameters`, `to_param`/`from_param_with_handler` (`param::DefaultParamHandler`, `DEFAULT_PARAM_HANDLER_SUPPORT.md`) |
| `ProgressLogger` base (`SignalToNoiseEstimator.h:32`): `startProgress` (:288), `setProgress` (:367), `endProgress` (:371), `setLogType` | `estimate_with_progress(..., &mut ProgressLogger)` and `estimate_peaks(..., Some(logger))`; the caller's `ProgressLogger::set_log_type`. The other entry points report nothing, as a source object whose log type is `NONE`. `NOISE_PROGRESS_LABEL` is the label. |
| `OPENMS_LOG_WARN` (:249, :379-382, :388-392) | `NoiseEstimates::log`, the lines without their line ends, as `RunOutput::log` in `feature_finder_picked` returns its lines |

### `SignalToNoiseEstimator.h` and `.cpp`

| Source (line) | Rust | Notes |
|---|---|---|
| class template `SignalToNoiseEstimator<Container>` (:30-33), `DefaultParamHandler` and `ProgressLogger` bases | trait `SignalToNoiseEstimator` | Implemented by `SignalToNoiseEstimatorMedian` and `SignalToNoiseEstimatorMeanIterative`. The parameter and progress surface is per type (above). |
| `typedef PeakIterator`, `typedef PeakType` (:39-40) | as rows 5 and 6 above | |
| constructor, copy constructor, `operator=`, destructor (:46-72) | `Default`, `Clone`, `Drop` of the implementors | |
| `virtual void init(const Container&)` (:75) | `SignalToNoiseEstimator::compute_stn`, `compute_stn_spectrum`, `compute_stn_chromatogram` | Returns the estimation instead of storing it. |
| `virtual double getSignalToNoise(Size index)` (:83) | `SignalToNoiseEstimator::signal_to_noise(&estimates)` | A slice. The source checks the index with `OPENMS_POSTCONDITION` only, which Release compiles out. The `@note` ("a warning to stderr if more than 20% ... sparse") is the median estimator's `:377` warning, returned in `NoiseEstimates::log`. |
| pure virtual `computeSTN_` (:96) | `SignalToNoiseEstimator::compute_stn` (required method) | `@exception InvalidValue` becomes `Error::InvalidValue`. |
| protected `struct GaussianEstimate { mean, variance }` (:105-109) | `GaussianEstimate { mean, variance }` | |
| protected `estimate_(first, last)` (:113-141) | `GaussianEstimate::of` (and crate-private `of_indexed`) | Population variance, sums in input order; an empty range gives the default NaN in both fields. |
| protected `stn_estimates_` (:146) | see State | |
| `estimateNoiseFromRandomScans(exp, ms_level, n_scans = 10, percentile = 80)` (:151; `.cpp:21-53`) | `estimate_noise_from_random_scans(exp, ms_level, n_scans, percentile, seed)`, `RandomScanNoise` | The seed is explicit (determinism contract). |
| `std::default_random_engine`, `std::uniform_real_distribution<double>` (`.cpp:34-35`) | `MinstdRand0`, `MinstdRand0::uniform01` | libstdc++'s `minstd_rand0` and `generate_canonical<double, 53>`, bit for bit. |
| `std::nth_element` (`.cpp:49`) | crate-private `noise_estimation::libstdcxx::nth_element` | The Release toolchain's algorithm, line by line. |

### `SignalToNoiseEstimatorMedian.cpp`

The file only defines the namespace-scope object
`SignalToNoiseEstimatorMedian<> default_sn_median2`, which no header declares;
it instantiates the `MSSpectrum` template inside libOpenMS and is not API. Not
ported; nothing can refer to it.

### Native additions

`NoiseEstimates::noise` and `max_intensity`; `max_points`, `max_bins`,
`max_work`; `RandomScanNoise::max_work` and `RANDOM_SCAN_NOISE_MAX_WORK`;
`NoiseCompatibility`, `BinIndexConversion`, `NoisePoint`;
`SIGNAL_TO_NOISE_ESTIMATOR_MEDIAN_NAME`, `NOISE_PROGRESS_LABEL`.

## Parameter contract

The nine defaults are `PEAK_PICKING_SUPPORT.md`'s `SignalToNoise:` list and
equal the stored `noise_defaults.ini`. Three bounds are open in the source and
therefore accepted under the source profile: `win_len` has no upper bound, so
`+inf` passes, and `NaN < 1` is false, so NaN passes too (`Param.cpp:146`,
the unset maximum being skipped at the same line); `noise_for_empty_window` has
no bound at all; `auto_max_stdev_factor` NaN passes `0..999` for the same
reason. An infinite `win_len` makes every window the whole container.

## Preserved source conventions

- **Histogram range.** `AUTOMAXBYSTDEV` is `sqrt(variance) * factor + mean`, as
  the Release build computes it, with both sums of `estimate_` in input order and
  the variance divided by `n`. `MANUAL` throws `Exception::InvalidValue` with the
  source text when `max_intensity <= 0`, only when estimation runs (:236-244).
  `AUTOMAXBYPERCENT` is below.
- **Negative range.** A range below zero writes the ungated warning
  `SignalToNoiseEstimatorMedian: the max_intensity_ value should be positive! <v>`
  (:249, `<v>` in `ostream` default format) and returns before the progress
  report starts, with every ratio `0.0` (the zero-initialised
  `stn_estimates_`), both percentages `0.0` and `noise` `+inf` (native field).
- **Windows.** The window of a point holds the points within `win_len / 2`
  (`0.5 * win_len` in the Release build) on either side, both ends inclusive.
  The two ends move as the source's iterators, which never lets the left end
  pass the right one: see [the window invariant](#the-window-invariant).
- **Bins.** The width is `max(1.0, max_intensity / bin_count)` (`maxsd`, so a
  NaN range gives width `1`). The bin of an intensity is
  `max(min((int)(I / width), bin_count - 1), 0)`, converted before clamping
  (CPP-257, below).
- **Median.** The median bin is the first whose cumulative count reaches
  `(count + 1) / 2`, never past the last bin, and the noise is
  `median * width + (half - before) / in_bin * width`, floored at `1.0`
  (`maxsd`, so NaN becomes `1.0`). A window with fewer than
  `min_required_elements` points is sparse and uses `noise_for_empty_window`.
  The ratio is `intensity / noise` with the intensity as the first operand.
- **Percentages.** `count * 100 / n`, not `count * (100 / n)`; the two differ
  for `3 * 100 / 7` and `5 * 100 / 6` (`pct_sparse_3of7`,
  `pct_rightmost_5of6`).
- **Warnings.** More than 20 % sparse windows and more than 1 % rightmost
  medians each add a line (:377-393), only when `write_log_messages` is set; the
  percentage is formatted as `std::ostream` does by default (`%g`, precision
  6). A NaN percentage compares false and writes nothing. The log stream's
  repeat suppression (`<line> occurred N times`) belongs to
  `concept::log_stream` and applies when the lines are written there.
- **Progress.** `startProgress(0, n, "noise estimation of data")` runs after the
  range is known and the early return was not taken, `setProgress(k)` after the
  `k`-th window and `endProgress()` after the loop, for an empty input too.
- **Release arithmetic.** Where IEEE-754 leaves a choice to the platform, the
  port makes the Release build's: a NaN operand of `addsd`/`subsd`/`mulsd`/
  `divsd` is returned quieted, the first operand's first; an invalid operation
  returns the negative default NaN `0xfff8000000000000` (`0xffc00000` in
  `f32`); `cvtss2sd` keeps a NaN's sign and payload. These are the private
  `noise_estimation::x86` functions, each use of which names the instruction.
- **Serial.** Neither header has OpenMP; the port is serial.

### AUTOMAXBYPERCENT on its defined domain

The source branch (:191-233), with `n = c.size()` `float` intensities `I_j`:

1. `:208-209`. `std::max_element` with the comparator `a > b` replaces its
   candidate whenever the candidate is greater than the next element, so it
   returns the **first minimum** `m` under `<` (the Release build emits
   `minss`, with the same NaN behaviour: a NaN candidate is never replaced, a
   NaN element never chosen). For an empty container it returns `end()`, and
   `:209` dereferences it: undefined (condition **D0**: `n >= 1`).
2. `:211`. `bin_size = m / 100`, a `float` division (`divss`), widened:
   `b = (double)(m / 100.0f)`.
3. `:216`. For every point, `q_j = (double)(I_j - 1.0f) / b` and
   `++histogram_auto[(int) q_j]` on a 100-element vector, unchecked. The write
   is in bounds exactly when the truncation of `q_j` lies in `[0, 99]`, i.e.
   `-1 < q_j < 100`, which excludes NaN. Outside, the write is out of bounds;
   where the truncation does not even fit `int`, the conversion (32-bit
   `cvttsd2si`) yields `INT_MIN`, which is out of bounds too. Condition **D1**:
   `-1 < q_j < 100` for every `j`.
4. `:220`. `t = (int)(p * n / 100)` with the integer `p` in `0..=100`, so
   `0 <= t <= n`: always defined.
5. `:221-230`. `i = -1`; while `run != end` and `seen < t`: `++i`,
   `seen += h[i]`, `++run`. Under D1 the 100 bins hold all `n >= t` points, so
   once `i = 99` has been read, `seen = n >= t` and the loop stops: it never
   reads `h[100]`. It also stops after `n` passes (one `run` step per pass),
   so `i <= min(99, n - 1)`, whatever the input order.
6. `:232`. `max_intensity = (i + 0.5) * bin_size`; `t = 0` leaves `i = -1` and
   a negative range, which takes the ungated early return.

So the branch is defined **exactly** on D0 and D1. What D1 means: if `m <= 0`,
either `b = 0` (a zero minimum, or a subnormal one that `m / 100` rounds to
zero), where `q` is infinite or NaN, or `b < 0`, where the minimum's own
quotient is `100 - 100 / m > 100`; so D1 implies `m > 0` and `b > 0`. Then the
minimum's quotient `100 (m - 1) / m` exceeds `-1` iff `m > 100 / 101`, and any
other point's quotient stays below `100` iff `I_j - 1 < m`, i.e. roughly
`I_j < m + 1`; NaN and `+inf` intensities are outside. The exact edges are the
`float` and `double` arithmetic above, which the port performs: the smallest
`f32` minimum of a constant input inside is `0.9900990128517151`
(`pctl_lower_edge`, the `f32` below it is outside), and for `m = 3` the largest
`f32` inside is `3.999999761581421` (`pctl_upper_edge`; `4` has quotient
`100.0000002`).

This refines CPP-256, which names the causes: the reversed comparator, a zero
bin size for a spectrum with a zero intensity, a negative index for intensities
below one, and the division by zero converted to `int`. Every one of those is
outside D0 or D1. The port computes the branch exactly on the domain
(`percentile_upper_end`) and returns `Error::Unsupported` elsewhere, naming
`:209` or `:216` and the point, in both profiles. Real spectra fall outside
almost always (the orbitrap spectrum does): the domain is a band of width one
above a minimum near or above one.

### CPP-257: the conversion before the clamp

`(int)(I / width)` at `:297` and `:308` is converted before `std::min<int>` and
`std::max` clamp it, so a quotient outside the `int` range, or a NaN, is
undefined. Both `libOpenMS.so` instantiations (spectrum and chromatogram) emit
the 32-bit `cvttsd2si %xmm1,%eax`, then `cmp %r12d,%eax; cmovg` and
`test %eax,%eax; cmovs` (`libOpenMS.so` at `0x186c5d0` and `0x186c64d` for
`MSSpectrum`). `cvttsd2si` returns the integer indefinite `INT_MIN` for NaN and
for every value whose truncation does not fit, which the clamp sends to **bin
0**. Measured on the Release build, twice per case, identically:
`cpp257_manual`, `cpp257_manual_chrom`, `cpp257_boundary` (`2^31 - 128`
still converts, `2^31` does not), `cpp257_stdev_bins`, `cpp257_neg_factor0`,
`cpp257_neg_default`, `cpp257_pick_noise`, and through `PeakPickerHiRes::pick`
(libOpenMS's own instantiation) `cpp257_pick_manual`, where the clamp-first
answer picks nothing and the Release build picks one peak.

The reach is wider than CPP-257 states: besides a manual `max_intensity` of
`1` with an intensity above `2^31`, the default standard-deviation mode reaches
it with `auto_max_stdev_factor = 0`, `bin_count = 1,000,000` and 10,000 points
with one at `f32(1e12)` (`cpp257_stdev_bins`: range `1e8`, width `100`,
quotient `1e10`), and negative intensities, which the source profile accepts,
pull the range down against the largest intensity (`cpp257_neg_factor0`,
`cpp257_neg_default`). Their own quotients can also fall below `INT_MIN`,
where both conversions give bin 0; only a quotient above `INT_MAX` or a NaN
separates them.

`BinIndexConversion::X86_64Release` (the source profile) reproduces the
measured outcome with `x86::cvttsd2si32`; `ClampBeforeTruncation` (the native
profile) clamps in `f64` first and follows the parameter description ("All
intensities EQUAL/ABOVE 'max_intensity' will be added to the LAST histogram
bin"). No case of the P1 oracle exercises this path at all:
`no_executed_case_bins_a_quotient_outside_the_int_range` checks every record of
every input of the 31 cases that run the estimator, and the P1 differential
runs every case with the clamp-first conversion as well and gets the same
bits.

### The empty-median-bin fallback is unreachable

`:350-353` falls back to the bin centre when the median bin is empty; the
source comment says this happens "if the rightmost bin was hit while empty". It
never happens. Every histogram count is the number of window points in that
bin, so the counts sum to `elements_in_window = count >= 1`, and
`half = (count + 1) / 2 <= count`. The walk (:325-329) enters a bin only while
the cumulative count before it is below `half`, and stops at the first bin
where the cumulative count reaches `half`, or at the last bin. If it stops by
reaching `half`, the stopping bin added at least one point. If it stops at the
last bin without reaching `half`, the cumulative count over all bins, `count`,
would be below `half <= count`, a contradiction. So the stopping bin always
holds at least one point, also for unsorted input and every bin conversion
(the add and the remove of a point use the same bin). The port keeps the branch
for fidelity with a comment saying so; no evidence of any tier can exist for
it. The earlier Rust comment repeated the source's wrong condition.

### The window invariant

The source's left end has no `left < right` guard (:295). It never needs one
on the parameter domain (`win_len >= 1`, `+inf` or NaN, so `h = win_len / 2` is
at least `0.5`, `+inf` or NaN):

- For a centre `c`, `pos(c) < pos(c) - h` is false for every position, so the
  left end stops at the centre at the latest.
- The right loop adds the centre when it reaches it, since
  `pos(c) <= pos(c) + h`, except when `pos(c)` is NaN, `pos(c) = -inf` with
  `h = +inf` (`-inf + inf` is NaN), or `h` is NaN. In those cases the right end
  stays at the centre, and at the following centres the left loop's condition
  on that point, `NaN < x`, `-inf < pos - inf` or `x < NaN`, is false again.

So the left end only removes points the right end added, `left <= right`
always holds, and the port's `left < right` guard never changes an iteration.
The oracle's `nf_pos_*` and `win_nan` cases pin the non-finite shapes.

### estimateNoiseFromRandomScans

Reproduced with its defined quirks (`RandomScanNoise::estimate`):

- The candidates are the non-empty spectra of `ms_level`; without one the
  result is `0.0` and `time` is not called (:32).
- The engine is `minstd_rand0` (`16807`, modulus `2^31 - 1`), seeded with the
  seed modulo the modulus, zero replaced by one. A draw is
  `generate_canonical<double, 53>`: two outputs `a`, `b` and
  `((b - 1) * r + (a - 1)) / r2` with `r = 2^31 - 2` and `r2` the double
  `0x43cfffffff000000`, clamped to `nextafter(1, 0)`; the constants are the
  `.rodata` values `libOpenMS.so` reads.
- The drawn index is `(UInt)(u * (candidates - 1))` (:42): 64-bit `cvttsd2si`,
  low 32 bits kept. The last candidate position is never drawn, and a single
  candidate always draws index `0`.
- :44 reads `exp[scan]`, not `exp[spec_indices[scan]]`, so the drawn spectrum
  may have any MS level (and may be empty).
- The position is `(Size)(size * percentile / 100.0)` (:48), converted as GCC
  converts a double to `unsigned long`: `comisd 2^63`; below it (and for NaN)
  a signed `cvttsd2si`, otherwise `subsd 2^63; cvttsd2si; btc 63`
  (`x86::f64_to_u64`). A value in `(-1, 0)` truncates to `0` (defined); a
  value of `2^64` or more, and `+inf`, gives `0` (the conversion is undefined
  there, and this is the Release outcome, `rnd_p_1e21`, `rnd_p_1e30`,
  `rnd_p_inf`); NaN gives `2^63`.
- `std::nth_element` (:49) runs the Release toolchain's algorithm
  (`libstdcxx::nth_element`, a line-by-line port of `bits/stl_algo.h` and
  `bits/stl_heap.h` of conda-forge GCC 14.4.0, sha256 `0598c5b1…` and
  `f18f83b2…`), so NaN intensities, for which `<` is not a strict weak ordering
  and the standard gives no guarantee, land where that build puts them. The
  algorithm stays in bounds for every irreflexive and asymmetric comparison,
  which `<` on `float` is with NaN; the argument is written at the module.
- The sum is `addss tmp[idx], noise` in draw order and the result
  `noise / (float) n_scans` (`divss`); `n_scans = 0` returns the negative
  default NaN.

## Native differences

The native safety profile differs from the source profile in exactly these
points (`NATIVE_REFUSALS`, `NATIVE_EMPTY` and `NATIVE_CLAMP` in
`tests/signal_to_noise.rs` list the oracle cases each one moves):

1. **Inputs.** Non-finite positions or intensities, negative intensities,
   duplicate positions (`Error::InvalidValue`) and decreasing positions
   (`Error::UnsortedData`) are refused. The source documents none of them as a
   precondition. `allow_negative_intensities`, `allow_duplicate_positions`,
   `allow_unsorted_positions` and `NoiseCompatibility::source_value_domain`
   lift them.
2. **Parameter values the source accepts.** `win_len` must be finite (the
   source accepts `+inf` and NaN), `noise_for_empty_window` finite and
   positive (the source accepts any value, including `0`, negative values,
   infinities and NaN, and divides by it), `auto_max_stdev_factor` not NaN.
3. **Non-finite results.** A non-finite histogram range (except for an empty
   input), window bound, noise or ratio is an error; the source stores what IEEE
   arithmetic gives, e.g. `+inf` for a positive intensity over a tiny positive
   `noise_for_empty_window` (`nfew_tiny`, `nfew_tiny_pos`).
4. **Empty input.** Zero for both percentages and, in the standard-deviation
   range, for `max_intensity`; the source divides zero by zero (:373-374 in
   every mode; `SignalToNoiseEstimator.h:127` in mode 0) and returns the
   default NaN, which `nan_for_empty_input` reproduces (`empty_stdev`,
   `empty_manual`, `empty_stdev_chrom`, `progress_empty`). Zeros keep a
   native caller's arithmetic finite.
5. **Bin conversion.** Clamp first (CPP-257, above).
6. **Resource ceilings.** `max_points` (1,000,000), `max_bins` (1,000,000),
   `max_work` (50,000,000 histogram updates and median-walk steps),
   `RandomScanNoise::max_work` (100,000,000 copied intensities plus draws).
   They apply in both profiles; the source has none.
7. **Explicit seed.** `estimateNoiseFromRandomScans` seeds with
   `time(nullptr)`; the port takes the seed as an argument.

Both profiles also differ from the source object in shape, not in results:

- **State: `stn_estimates_` and the percentages.** The source object stores the
  last estimation in `stn_estimates_`, `sparse_window_percent_` and
  `histogram_oob_percent_`; the port returns them as a `NoiseEstimates` value
  and holds no state, so nothing can be stale. Two source defects that this
  removes: neither constructor initialises the two percentage members (the
  default constructor's `defaultsToParam_` runs `updateMembers_`, which does
  not set them, and the copy constructor does not copy them), so calling a
  getter before `init` reads an indeterminate `double`; and `operator=` keeps
  the target's old percentages while it clears the estimates (CPP candidate,
  as CPP-289).
- **Logging.** The warnings are returned in `NoiseEstimates::log` instead of
  being written to `OPENMS_LOG_WARN`.
- **Progress.** Opt-in through `estimate_with_progress`; a native refusal after
  the report started ends the report before returning, so the logger's nesting
  stays balanced (the source cannot fail there).

## Refusals in the source profile

Each is exactly as wide as the undefined behaviour it avoids.

| Source | Undefined behaviour | Refused domain |
|---|---|---|
| `SignalToNoiseEstimatorMedian.h:209` | `*end()` of an empty container | `auto_mode 1` with no points |
| `:216` | write outside the 100-element pre-histogram | `auto_mode 1` with any point's `q` outside `(-1, 100)` (the conversion's `INT_MIN` included) |
| `:365` (`++window_count`), `SignalToNoiseEstimator.h:123` (`++size`) | signed `int` overflow | more than `i32::MAX` points (the native `max_points` ceiling is lower by default) |
| `:324` (`elements_in_window + 1`) | signed `int` overflow | a non-sparse window of exactly `i32::MAX` points |
| `SignalToNoiseEstimator.cpp:49` | `tmp.begin() + idx` past `end()` | a drawn scan whose position exceeds its size: a percentile above `100`, a product of `-1` or below, a NaN product |
| `SignalToNoiseEstimator.cpp:50` | `tmp[idx]` one past the end | position equal to the size: percentile `100`, an empty drawn scan |

Out-of-domain probes on the Release build (`../oracle/sne-completion/probes/`,
each run three times) show there is no answer to reproduce: seven of the nine
`auto_mode 1` probes end with `SIGSEGV` (139) or `SIGABRT` (134) every time;
the other two (a write exactly one bin past the end, and a negative minimum)
return, as does `percentile 100` and `150`, with values read from neighbouring
memory (`0x6e`, `0x21`); an empty drawn scan ends with `SIGSEGV`. One probe is
deterministic in a way the table's rule still refuses: a **NaN percentile**
returned the scan's minimum in all three runs, because the index `2^63` times
four wraps to the first element in `lea (%r12,%rcx,4)`. The C++ pointer
arithmetic at `:49` is undefined all the same, and the lead's rule for this
wave limits emulation to float-to-integer conversions, so the port refuses it
and records the measurement here (the same wrap would read element `k` for any
position `2^62 j + k`, which only a product of at least `2^62` can produce).

## Undefined behaviour reproduced (x86-64 Release)

| Site | Measured outcome | Instruction | Emulation | Cases, repetitions |
|---|---|---|---|---|
| `SignalToNoiseEstimatorMedian.h:297`, `:308` | bin 0 for NaN and out-of-range quotients | 32-bit `cvttsd2si`, then `cmovg`/`cmovs` clamp | `x86::cvttsd2si32` via `BinIndexConversion::X86_64Release` | 8 `cpp257_*` cases, each run twice |
| `SignalToNoiseEstimator.cpp:48` | position `0` for products of `2^64` or more and `+inf` | `comisd 2^63; subsd 2^63; cvttsd2si; btc 63` | `x86::f64_to_u64` | `rnd_p_1e21`, `rnd_p_1e30`, `rnd_p_inf`, each twice |
| `SignalToNoiseEstimator.cpp:42` | low 32 bits of the truncation (only past `2^32` candidates, unreachable in memory) | 64-bit `cvttsd2si; mov %esi,%esi` | `x86::cvttsd2si64(..) as u32` | the 52 random-scan cases stay in range |
| `SignalToNoiseEstimator.cpp:49` (not a conversion) | the Release library's permutation for NaN input | `std::nth_element` of GCC 14.4.0 | `libstdcxx::nth_element` | `nth_vectors` (390 vectors, 30 through `__heap_select`), 19 `rnd_nan_*` cases, each twice |

## Class-test accounting

`SignalToNoiseEstimatorMedian_test.cpp` (5 sections): the constructor, copy
constructor, assignment and destructor are `Default`, `Clone` and `Drop`
(`NOT_TESTABLE` in the source; the copy sections call `init` on an empty
spectrum with the defaults, which the oracle's `empty_stdev` executes);
`[EXTRA] init` with `win_len 40`, `noise_for_empty_window 2`,
`min_required_elements 10` against `SignalToNoiseEstimatorMedian_test.out` is
`class_test_noise_estimator_section` (`tests/peak_picking_experiment.rs`,
`TEST_REAL_SIMILAR`), and bit for bit the P1 oracle's `class_noise_init` and
this oracle's `dta_class`.

`SignalToNoiseEstimator_test.cpp` (6 sections): the four special members and
`init` of a test subclass whose `computeSTN_` does nothing, and
`getSignalToNoise` (`NOT_TESTABLE`):
`class_test_base_sections_through_a_trivial_estimator`
(`tests/signal_to_noise.rs`).

No class test calls `estimateNoiseFromRandomScans`.

## Evidence

- **Tier 1, Linux x86-64 Release, direct.** `tests/signal_to_noise.rs`
  compares the 143 cases of `tests/data/signal_to_noise/cases.tsv` (over
  `synthetic.tsv`, the class-test DTA and two P1 mzML inputs) with the
  unmodified Release build `openms4-release-bc9cc12-c19e494-174b576` on
  `ibminode06`, each case in a fresh process, run twice byte-identically
  (driver, generators, disassembly, probes and hashes in
  `../oracle/sne-completion/`; `tests/data/signal_to_noise_provenance.json`).
  All values are IEEE-754 bit patterns.
  - 89 estimator cases, 31,322 ratios: `AUTOMAXBYSTDEV` 48, `MANUAL` 22 and
    `AUTOMAXBYPERCENT` 19 (the parameter's value, counted per case), with
    every ratio, `max_intensity_`, both percentages and the count of
    estimates.
  - Percentages strictly between 0 and 100 in 8 cases, including the three
    where `count * 100 / n` and `count * (100 / n)` differ
    (`pct_sparse_3of7`, `pct_rightmost_5of6`, `pct_both_3of14`).
  - `SignalToNoiseEstimatorMedian<MSChromatogram>` directly in 8 cases
    (`chrom_*`, `cpp257_manual_chrom`, `pctl_const_chrom`,
    `empty_stdev_chrom`, `progress_chrom`).
  - The negative-range early return in 7 cases, and the three warnings and
    their gating byte for byte against the C++ standard error (`warn_sparse`,
    `warn_clipped`, `warn_negative_range`, the `_nolog` twins,
    `pct_both_3of14`).
  - The progress report (`progress_*`, 7 cases) byte for byte against the C++
    standard output with the two timings masked; `time()` is interposed by the
    driver, constant or advancing by one per call, which makes the
    once-per-second throttle deterministic.
  - The two exceptions, with their text (`manual_invalid`, `manual_zero`).
  - `estimateNoiseFromRandomScans` in 47 cases with the seed set through the
    interposed `time()` (one call per estimation that has candidates, none
    otherwise), the engine's raw outputs and uniform draws in 5, and the whole
    `nth_element` permutation of 390 vectors in 1
    (`permutations_match_the_release_toolchain`).
  - The tool path: `cpp257_pick_manual` runs `PeakPickerHiRes::pick`, i.e.
    libOpenMS's own instantiation.
- **Which instantiation.** The estimator is a header-only template, so the
  driver instantiates it itself, with the same compiler and libOpenMS's own
  compile flags (`-O3 -DNDEBUG -std=gnu++23 -fPIC -fvisibility=hidden
  -fvisibility-inlines-hidden -mssse3 -ffp-contract=off -fopenmp`, the entry
  for `SignalToNoiseEstimator.cpp` in the Release build's
  `build/core/compile_commands.json`). Its 88 floating-point instructions
  equal libOpenMS's in order and operands for `MSChromatogram`, and for
  `MSSpectrum` differ in one memory displacement only (`cvtss2sd -0x8(REG)`
  against `0x8(REG)`, a different pointer increment), with the same four
  32-bit `cvttsd2si`. The tool-path case exercises libOpenMS's copy directly.
- **Tier 1, the P1 oracle, direct and indirect.** `tests/peak_picking_experiment.rs`
  compares 93 cases with the **arm64 Debug** product SDK (core `4fdec46`, whose
  diff to `bc9cc12` is empty for these sources). Of them, 11 cases carry
  estimator ratios and percentages directly (the 12 `noise` cases less
  `extra_noise_manual_invalid`), and their percentages are only ever `0` or
  `100`; 19 picker cases run the estimator indirectly (`signal_to_noise > 0`),
  among them 6 chromatogram records (5 in `extra_topp2_parameters`, 1 in
  `extra_topp2_chromatogram0_check`); the remaining 62 cases never run it,
  because `signal_to_noise` defaults to `0` (`PeakPickerHiRes.cpp:31`).
  The unchanged P1 driver (sha256 `2d06db2f…`) re-run against the Linux
  x86-64 Release build printed byte-identical records (sha256 `9eb8f249…`)
  and parameter files; its standard error is identical except for six lines
  that only the Debug build prints, five mzML loader diagnostics and one
  `Update ranges was called but ranges were already up-to-date`
  (`../oracle/sne-completion/p1/`).
  So the P1 fixture is also the Release build's output.
- **Tier 3.** The class-test sections above.
- **Tier 4.** Native refusals and their exact case lists; the `x86`
  emulations' edge values; the `libstdcxx` port's order-statistic and
  multiset invariants on 20,000 random vectors with and without NaN; progress
  balance on a refusal; parameter round trips, the handler's `param_`, and the
  gated and ungated warnings; the base trait on both implementors.

## Open points

- `PeakPickerChromatogram` (`src/processing/chromatogram.rs:181`, another
  header) calls the strict `estimate`; the source accepts the duplicate and
  negative samples a smoothed chromatogram can have.
- The NaN-percentile wrap above is refused under this wave's rule; emulating it
  would need the rule to cover pointer arithmetic.
