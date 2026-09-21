# Iterative peak picking

`processing::iterative::PeakPickerIterative` ports the algorithm in
`PROCESSING/CENTROIDING/PeakPickerIterative.h` at OpenMS4-core revision
`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The `.cpp` contains no algorithm;
the source class tests exercise construction and leave numerical tests as TODOs.
The native numerical tests therefore distinguish analytical examples and
independent source-derived calculations from preexisting C++ numerical goldens.

## API and defaults

```rust
use openms::kernel::{MSSpectrum, Peak1D};
use openms::processing::iterative::PeakPickerIterative;

let input = MSSpectrum {
    peaks: [0.0, 1.0, 4.0, 10.0, 4.0, 1.0, 0.0]
        .into_iter().enumerate()
        .map(|(i, y)| Peak1D::new(100.0 + i as f64, y)).collect(),
    ..Default::default()
};
let picker = PeakPickerIterative { signal_to_noise: 0.0, ..Default::default() };
let result = picker.pick_spectrum(&input)?;
assert_eq!(result.picked.spectrum.peaks, [Peak1D::new(103.0, 20.0)]);
# Ok::<(), openms::Error>(())
```

The defaults are signal-to-noise 1, peak half width 0, spacing multiplier 1.5,
five iterations, internal width checking disabled, all MS levels selected, and
generated float arrays retained. Refinement uses a median-noise estimator with
30 bins and a 20-unit window. The typed `noise_estimator` permits the existing
estimator's other settings; its defaults otherwise match OpenMS.

HiRes seed picking receives the same signal-to-noise and spacing settings, but
retains its **independent default noise settings**, including a 200-unit window
and 30 bins. Altering refinement noise does not silently change seed detection.
Signal-to-noise zero disables both estimators. The seed stage retains the
already documented native HiRes checks and numerical policies.

`pick_spectrum` returns `IterativePickingResult` with a `PickedSpectrum` and
aligned `IterativePeakRegion` entries. Each region reports the seed ordinal,
initial raw center index, final recentered raw index, and inclusive left/right
raw indices. `picked.boundaries` stores those final boundaries at full `f64`
precision. `picked.omitted_arrays` names all input profile arrays without an
aggregation rule. Record metadata is preserved, including retention time,
native ID and identification records; representation becomes centroid.

Three output float arrays retain the source names and order:

| Array | Contents |
|---|---|
| `IntegratedIntensity` | Inclusive sum of raw intensities, also used as output intensity |
| `leftWidth` | Original left boundary m/z rounded to `f32` |
| `rightWidth` | Original right boundary m/z rounded to `f32` |

These are sampled intensity sums, not trapezoidal areas. All input float,
integer and string arrays are omitted and reported, including ion mobility.
Seed detection uses only raw m/z and intensity; it does not compute unused
mobility aggregates. Full precision boundary values and original indices remain
available even when float rounding changes their numerical values.

`pick_experiment` returns `IterativeExperimentResult`. It preserves
chromatograms and experiment metadata. `ms1_only` copies excluded spectra
unchanged; their `spectrum_regions` entry is `None`. `clear_meta_data` clears
the three generated float arrays only for spectra picked through this experiment
method, leaving exact regions in the result. These two options do not affect
direct `pick_spectrum` calls. `SpectrumFilter` provides atomic spectrum and
experiment mutation wrappers.

The source's `ProgressLogger` base, as `pickExperiment` uses it, is
`pick_experiment_with_progress(input, &mut ProgressLogger)`: the caller passes
the logger for the call and selects its type on it. It makes the source's calls
(`PeakPickerIterative.h:384-398`): one section over the spectrum count labelled
`picking peaks` (`ITERATIVE_PICKING_PROGRESS_LABEL`), advanced after every
spectrum, picked or copied, with the value *before* the source's
post-increment, so the values run from `0` to `n - 1` where `PeakPickerHiRes`
reports `1` to `n`. Chromatograms, which the port keeps, are not counted,
because the source's range does not count them. The picked result is the one
`pick_experiment` returns. `tests/progress_consumers.rs` replays the Release
build's calls and command output for three spectra and for an empty experiment
(tier 1). Validation and the point-limit preflight run before the section, so
an input refused there prints nothing; a failure inside the section still ends
it, which the source, whose loop cannot fail there, never needs to do.
`pick_experiment` itself reports nothing, as a source object of type `NONE`.

## Preserved refinement behavior

Seed association scans raw samples in increasing m/z and assigns each seed to
the first subsequent raw point **strictly greater** than its coordinate. It
advances at most one seed per raw sample. It does not choose the nearest raw
point or perform an independent upper-bound search for every seed.

Candidates are prioritized by the original HiRes seed intensity, not their
later integrated intensity. Each iteration includes the current raw center and
its immediate neighbors unconditionally. Optional internal width checking
discards a candidate when either immediate spacing is greater than `peak_width`;
equality passes. Enabling it while keeping zero width removes all candidates
on a valid profile with distinct coordinates, as in the source.

Further extension requires a spacing strictly below `spacing_difference`
times the smaller initial neighboring spacing. That initial spacing remains
fixed throughout the extension. Intensities must strictly decrease, unless
the sample is strictly inside the configured half width. Active noise checks
stop when the ratio is below the threshold; equality passes. The central three
samples are not subjected to this refinement noise gate.

Integration accumulates samples from left to right in `f64`, matching the
source's ordered map traversal. The intensity-weighted centroid is rounded to
`f32` for candidate storage, overlap suppression and output, then promoted to
the kernel's `f64` m/z field. Raw recentering instead uses the unrounded
weighted mean. It searches the left side first, changes center only on a
strict improvement, and excludes boundary samples and raw index zero.

The source rightward search uses `i - m > 0` while accessing `i + m`. This
asymmetry is retained, with an additional upper-index guard. Consequently it
can stop before reaching the globally nearest raw point. Subsequent iterations
use this constrained center. No symmetric nearest-point replacement is made.

After all iterations, each surviving candidate suppresses lower-priority
candidates whose rounded centroids lie inclusively inside its full-precision
boundaries. Original seed priority remains unchanged. Surviving peaks, regions
and arrays are sorted together by their rounded output m/z. Equal seed
priorities retain original seed order; this deterministic choice replaces
the C++ `std::sort` unspecified ordering of equivalent elements. Equal rounded
output coordinates retain that candidate-priority order.

## Which profile the noise estimate uses

Source `pick` builds a fresh `SignalToNoiseEstimatorMedian<MSSpectrum>`, sets
`win_len` and `bin_count` on it and calls `snt.init(input)` over the caller's
raw spectrum, but only when `signal_to_noise_` is positive
(`PROCESSING/CENTROIDING/PeakPickerIterative.h:314-321`). Its seed picker `pp`
is likewise a default-constructed `PeakPickerHiRes` with two parameters changed
(`:288-292`). Both are plain source objects, and the signal they read is the
caller's, so `PeakPickerIterative::compatibility` is handed to both, exactly as
`PeakPickerHiRes` uses its own.

`compatibility` lifts one refusal here:

| Flag | Source behaviour |
|---|---|
| `allow_negative_intensities` | No intensity check exists in the source. Negative samples reach `snt.init`, and a candidate whose support sums to a negative intensity divides by that sum (`PeakPickerIterative.h:228`) to a finite recentred m/z and is stored with that negative intensity (`:231-234`). Measured as `ppi_neg_sn0`, `ppi_negbase_sn0` and `ppi_allneg_sn0`, which reach the negative sum because `signal_to_noise_ = 0.0` turns the S/N gate off (`:92`, honoured at `:185`, `:205` and `:315`). A zero sum stays refused in both profiles: it divides to an infinite or NaN centroid, and no case pins the order `std::stable_sort` then produces. |

The `noise` sub-profile carries the estimator's own source behaviours, including
the `win_len` values the source's `setMinFloat("win_len", 1.0)` restriction lets
through — NaN, because `NaN < 1` is false, and `+inf`, because the unset upper
bound is skipped — which the native median-noise profile refuses. Values the
restriction rejects (`0.5`, `-inf`, `bin_count` of `1` or `2`) throw
`Exception::InvalidParameter` in the Release build before `init` runs, and are
refused in both profiles here.

Reading the peaks through `estimate_peaks` rather than copying two `f64` slices
also fixes the intensity widening: the source's `getIntensity()` feeds
`computeSTN_` through `cvtss2sd`, and `x86::widen` reproduces that for every NaN
payload, where a `f32`-to-`f64` cast in Rust is free to produce any NaN.

### Two flags this picker does not yet honour

`allow_duplicate_positions` and `allow_unsorted_positions` leave the
corresponding refusals in place in both profiles. The source does reach
`snt.init` with such spectra — measured, the Release build picks both without
complaint — but `pickRecenterPeaks_` then collects each peak's support in a
`std::map<double, double>` keyed by m/z: an equal m/z overwrites the stored
intensity and contributes once to the integrated intensity and to the weighted
centroid, the spacing bounds come from `begin()`/`rbegin()` rather than from the
first and last visited index, and every spacing is taken through `std::fabs`.
This port integrates over an index range, which agrees with all of that only for
strictly increasing positions. Accepting either flag before those semantics are
ported would return different peaks under a flag that claims source behaviour,
so both stay refused. Porting them is tracked as remaining scope below.

## Checked differences and resource bounds

- Profile m/z must be finite, nonnegative and strictly increasing; intensities
  must be finite, and nonnegative unless
  `PickingCompatibility::allow_negative_intensities` is set. Negative
  coordinates would collide with source `-1` deletion sentinels. Duplicate
  positions would be collapsed by the C++ map and are rejected instead. Zero
  intensities are accepted.
- Zero iterations, zero limits, negative or nonfinite width/noise/spacing
  settings, missing center neighbors and nonfinite arithmetic return errors. A
  nonpositive integrated intensity returns an error in the native profile; under
  `allow_negative_intensities` a negative sum divides as the source divides by
  it (`PeakPickerIterative.h:228`) and only a zero sum still fails, through the
  infinite or NaN centroid it produces. Values not representable as finite output `f32`
  values also return errors. A finite rounded output centroid can lie outside
  its exact integration boundaries; the source rounding behavior is retained.
- Empty and short profiles return empty centroid spectra with preserved record
  metadata and the three empty output arrays. C++ returns without assigning
  output for fewer than three points. Native experiment picking also preserves
  input chromatograms, rather than dropping them when rebuilding the map.
- Per-spectrum `max_points` defaults to 1,000,000 and is checked before input
  validation and allocations. Selected experiment spectra are preflighted before
  the experiment is cloned. Copied records and record metadata retain their
  existing sizes; this is not a global experiment-memory limit.
- `max_work` defaults to 50,000,000 per stage. It bounds HiRes candidate/support
  visits, each active noise histogram's setup plus estimator visits, and one
  combined budget for association, all refinement iterations, integration,
  recentering, sorting allowance and quadratic overlap comparisons. The sorting
  allowance is `2*c*ceil(log2(c))` for `c` candidates. Child estimator limits
  are respected when smaller. Histogram setup is charged before allocation.
  Huge iteration counts are rejected before iteration work when candidates
  exist; empty candidate sets do not execute an empty iteration loop.

With `n` samples, `c <= n` seeds and `r` iterations, refinement is worst-case
`O(r*c*n + c*c + c*log(c))`; ordinary narrow peaks visit fewer samples.
Working storage is `O(n + c + b)` for `b` active histogram bins, plus copied
record metadata and annotations. Existing HiRes and median-noise stage bounds
apply separately. All APIs calculate results before replacing caller data, so
a failure in a later spectrum does not leave an experiment partly processed.

## Verification

Focused public tests cover an analytical symmetric peak, raw intensity sums,
default sparse-noise rejection, internal width equality, original boundaries,
metadata and omitted arrays, experiment-only flags, unchanged chromatograms,
atomic failure, malformed profiles and resource limits. Independent private
tests inject explicit seeds into the natural refinement helper and compare
source-derived cases for assignment, strict thresholds, priorities, inclusive
suppression, right-search asymmetry and `f32` rounding. External reference tests
also exercise the public picker on pinned Orbitrap/FTMS profile inputs.
Fixture provenance records the source and independent extraction semantics;
these checks are not labeled as nonexistent upstream iterative output goldens.
