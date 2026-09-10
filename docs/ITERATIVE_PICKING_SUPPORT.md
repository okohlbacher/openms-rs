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

## Checked differences and resource bounds

- Profile m/z and intensities must be finite and nonnegative, with strictly
  increasing m/z. Negative coordinates would collide with source `-1` deletion
  sentinels. Duplicate positions would be collapsed by the C++ map and are
  rejected instead. Zero intensities are accepted.
- Zero iterations, zero limits, negative or nonfinite width/noise/spacing
  settings, missing center neighbors, zero integrated intensity and nonfinite
  arithmetic return errors. Values not representable as finite output `f32`
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
