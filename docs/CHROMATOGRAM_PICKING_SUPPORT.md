# Chromatogram peak picking

`processing::chromatogram::PeakPickerChromatogram` ports the native `legacy` and
`corrected` paths of pinned
[`PeakPickerChromatogram.cpp`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp).
It reuses the existing Gaussian/Savitzky–Golay filters, `PeakPickerHiRes`, and
histogram-median noise estimator. Crawdad is an external backend and is not
implemented or represented by either native method.

## API and defaults

```rust
use openms::kernel::MSChromatogram;
use openms::processing::chromatogram::PeakPickerChromatogram;

let input = MSChromatogram::default();
let result = PeakPickerChromatogram::default().pick_chromatogram(&input)?;
let peaks = &result.picked.chromatogram;
for region in &result.regions {
    let left_rt = input.peaks[region.left_index].rt;
    let right_rt = input.peaks[region.right_index].rt;
    // Both endpoints are original samples, included in the raw intensity sum.
    assert!(left_rt <= right_rt);
}
# Ok::<(), openms::Error>(())
```

`pick_chromatogram(&MSChromatogram)` returns `ChromatogramPickingResult` with:

- `picked: PickedChromatogram`: spline apices, source output arrays, the HiRes
  seed-support `boundaries`, and `omitted_arrays` identifying profile arrays
  without a centroid aggregation rule.
- `smoothed: MSChromatogram`: the complete original sampling and annotations,
  with smoothed intensities.
- `regions: Vec<ChromatogramPeakRegion>`: `apex_index`, `left_index`, and
  `right_index`, all referring to the original input. The boundary indices are
  inclusive and describe the final integration regions. They differ from the
  HiRes seed support stored in `picked.boundaries`.

`filter_chromatogram(&mut MSChromatogram)` commits only the picked chromatogram
after the entire operation succeeds. Use the returning method to retain the
smoothed trace, original indices, and omitted-array report.

| Option | Default and meaning |
| --- | --- |
| `method` | `ChromatogramPickingMethod::Corrected`; `Legacy` is also supported |
| `smoothing` | `ChromatogramSmoothing::Gaussian { width: 50.0 }`, in seconds |
| `peak_width` | `None`; `Some(positive_seconds)` enables forced extension |
| `signal_to_noise` | 1.0, for boundary extension |
| `noise_estimator` | Median estimator, 1000-second window and 30 bins |
| `seed_signal_to_noise` | 1.0, independently controlling HiRes seeds |
| `seed_noise_estimator` | Median estimator, 200-second window and 30 bins |
| `report_sn` | false; enables reporting when boundary S/N is disabled |
| `remove_overlapping_peaks` | false |
| `max_points` | 1,000,000 |
| `max_work` | 50,000,000 per stage, as detailed below |

The alternative smoothing option is
`ChromatogramSmoothing::SavitzkyGolay { frame_length: 15, polynomial_order: 3 }`,
matching the source's configured SG defaults. Even frame lengths are incremented.
SG fits sample indices rather than RT distances; the existing filter's numerical
limits apply. Gaussian width is eight standard deviations, and its existing
0.01-second kernel lookup spacing and endpoint conventions are preserved.

## Separate seed and boundary noise settings

The C++ constructor configures its internal HiRes picker once, with S/N 1.0 and
HiRes's default 200-second noise window. Later changes to the outer
`signal_to_noise`, `sn_win_len`, or `sn_bin_count` do not update that seed picker.
The native API preserves these initial values and exposes both configurations
explicitly. Changing `signal_to_noise` to zero therefore does not disable seed
filtering; change `seed_signal_to_noise` separately to do that.

Seed picking always uses the smoothed trace, with spacing constraints disabled
and FWHM reported in absolute seconds. The returned peak positions and intensities
are the spline apices of that smoothed signal. `Legacy` chooses integration bounds
on the original intensities; `Corrected` chooses them on the smoothed intensities.
Both methods integrate the original raw samples.

Boundary S/N is estimated on the same trace used to choose bounds. If
`signal_to_noise > 0`, apex S/N values are reported even when `report_sn` is false.
If the threshold is zero and `report_sn` is true, noise estimation is used only
for the report. If both are disabled, every output `SN` value is -1. The apex
report uses the closest original sample index, not the interpolated intensity.

## Boundary and overlap semantics

The source's closest-sample lookup walks forward from the previous apex, chooses
the right sample at an exact distance tie, and returns the end index for targets
at or beyond the final sample. The native caller checks the end sentinel and
requires neighbors on both sides of each seed. It does not index an invalid
endpoint or reproduce unsigned-index underflow.

Each region initially includes the immediate samples on both sides of the apex.
Tests begin at the second sample away. A candidate extends the region when its
intensity is strictly below its inward neighbor, or when its distance from the
spline apex is strictly below `peak_width`. In either case it must pass the
configured S/N threshold when that threshold is positive. Forced width does not
override the noise gate. A failed candidate is not included, while the immediate
neighbors remain included unconditionally.

Overlap adjustment compares adjacent regions and acts only when the first right
index is greater than the next left index. It walks down the intensity slopes
from both apices to find a valley. If the resulting bounds still cross, it
preserves the C++ sequential integer assignments:

```text
new_left  = floor((old_left + old_right) / 2)
new_right = floor((new_left + old_right) / 2)
```

The second expression uses the updated left value. These can leave overlapping
regions; this option does not promise a partition into disjoint intervals.
Even equal final bounds share their endpoint because integration is inclusive.
The source apices and their FWHM values are not replaced by the boundary search.
For mismatched legacy raw/smoothed shapes, the sequential adjustment can even
move a boundary past its unchanged closest-apex index.

## Five output arrays and precision

Arrays are stored as f32 in source order:

| Name | Meaning |
| --- | --- |
| `FWHM` | Width in seconds from the initial HiRes seed picker |
| `IntegratedIntensity` | Inclusive sum of original sample intensities, accumulated as f64 then stored as f32 |
| `leftWidth` | RT of the first included original sample, cast to f32 |
| `rightWidth` | RT of the last included original sample, cast to f32 |
| `SN` | Boundary-estimator S/N at the closest original apex sample, or -1 |

`IntegratedIntensity` is a sample sum. It does not multiply by elapsed time or use
trapezoidal integration. Irregularly sampled traces follow the same rule.
Use a separate integration method when a time-weighted area is required.

The left/right metadata values can round away from an input RT. For subsequent
integration, use `input.peaks[region.left_index].rt` and
`input.peaks[region.right_index].rt` to preserve the original full-precision
endpoints. All f32 conversions are checked; a finite f64 RT or sum that cannot
be represented as finite f32 returns an error rather than creating infinite
metadata.

The native kernel's chromatogram name, native ID, precursor and metadata survive
on both returned traces. All profile arrays stay aligned on `smoothed`; the
picked trace reports their names through the existing HiRes `omitted_arrays`.
No profile-array reduction is invented.

## Which profile the noise estimate uses

The C++ picker owns one `SignalToNoiseEstimatorMedian<MSChromatogram>` member
and configures it with `win_len`, `bin_count` and `write_log_messages` only
(`ANALYSIS/OPENSWATH/PeakPickerChromatogram.cpp:408-412`). It calls
`snt_.init(chromatogram)` at `:171` on the boundary signal, which is the
caller's chromatogram under `legacy` and the smoothed trace under `corrected`.
There is no native variant of that member: by the time it runs, the caller's
chromatogram has already passed `pickChromatogram`'s own checks, and under
`corrected` the signal is one this picker produced. The port therefore estimates
with `PickingCompatibility::source()` at that call site regardless of
`compatibility`, and `compatibility` governs only what the picker accepts from
its caller.

`compatibility` lifts exactly two refusals, both measured against the Release
build:

| Flag | Source behaviour |
|---|---|
| `allow_negative_intensities` | The source has no intensity check anywhere. A baseline-subtracted chromatogram under `legacy` puts negative samples straight into `snt_.init`, which bins them in the first histogram bin or, when the automatic range comes out negative, reports zero everywhere. |
| `allow_duplicate_positions` | `MSChromatogram::isSorted` accepts equal retention times, so they reach both the seed picker and `snt_.init`. |

`allow_unsorted_positions` is deliberately inert here. `pickChromatogram`
(`PeakPickerChromatogram.cpp:68-72`) throws `Exception::IllegalArgument`,
"Chromatogram must be sorted by position", so decreasing retention times stay
`Error::UnsortedData` in both profiles. `MSChromatogram::isSorted` decides that
with `prev.getRT() > next.getRT()` (`KERNEL/MSChromatogram.cpp`), which is why
equal retention times pass it, a NaN passes it (every comparison against a NaN
is false), and an infinity ahead of a finite sample does not — the measured
`IllegalArgument` for the `nonfinite` case comes from that last one.
Non-finite retention times and intensities are refused earlier here, by
`MSChromatogram::validate`, and no flag lifts that; `PeakPickerHiRes` refuses
them unconditionally too.

Two smoother facts decide how far a negative sample travels, and both were read
from the pinned source rather than assumed. `SavitzkyGolayFilter`
(`PROCESSING/SMOOTHING/SavitzkyGolayFilter.h:115`, `:135`, `:153`) writes
`std::max(0.0, help)`, so an SG-smoothed trace is never negative. `GaussFilter`
(`PROCESSING/SMOOTHING/GaussFilter.cpp:138`) does not clamp, but its kernel is
non-negative, so it only carries a negative sample through where the local
weighted average is itself negative. Under `corrected` the estimator therefore
usually sees a non-negative trace even from a negative input; under `legacy` it
sees the negatives directly.

The parameter domain follows the source's own `Param` restrictions, which
`DefaultParamHandler::setParameters` enforces through `Param::checkDefaults`.
`setMinFloat("win_len", 1.0)` rejects `0.5` and `-inf` with
`Exception::InvalidParameter` but passes NaN, because `NaN < 1` is false, and
`+inf`, because the unset upper bound is skipped; `setMinInt("bin_count", 3)`
rejects `1` and `2`. The port refuses and accepts exactly the same set, and a
NaN or infinite window picks rather than failing because the estimate uses the
source profile.

Using the source profile at that call site also selects the Linux x86-64
Release build's bin-index conversion, which is the one behavioural change the
default `compatibility` sees beyond the non-finite window. The source truncates
before it clamps — `std::max(std::min<int>((int)(I / bin_size), bin_count - 1),
0)` at `SignalToNoiseEstimatorMedian.h:297` and `:308` — so a quotient that
leaves `int` range becomes `INT_MIN` in that build and lands in bin `0`, where
the port's native conversion clamps it into the last bin. The quotient is
`I / std::max(1.0, max_intensity_ / bin_count_)` (`:258`), so producing an
out-of-range one needs `max_intensity` set by hand, which the source's picker
never does: `updateMembers_` sets only `win_len`, `bin_count` and
`write_log_messages` (`PeakPickerChromatogram.cpp:408-412`), leaving
`max_intensity` at `-1` and the range automatic at `mean + 3 sd`. With that
automatic range the bin width scales with the data — for one sample of `I`
among `n` near-zero ones the quotient is about `sqrt(n) * bin_count / 3`, so
exceeding `2^31` over 48 samples would take a `bin_count` near `9.3e8` and a
histogram of about eleven gigabytes. The Rust `noise_estimator` field exposes the whole estimator,
including `histogram_range`, so a caller can build the configuration the source
cannot; there the estimate is the Release build's own, measured as the
`ppc_bigmax` / `ppi_bigmax` cases in `tests/data/picker_consumers/snt_oracle.tsv`
(a hand-set upper end of 10 over 30 bins floors the bin width at 1.0, so a 3e9
sample's quotient leaves `int` range) and asserted in
`tests/picker_noise_consumers.rs`.

## Validation, limits and remaining scope

Coordinates must be finite and non-decreasing; intensities must be finite.
Decreasing retention times are refused in both profiles, as the source throws;
equal ones need the source profile.
Duplicate RT samples and negative intensities are refused by the native profile
and accepted by the source profile, as the table above records. Malformed
parallel arrays, invalid active noise parameters, nonpositive forced widths,
unsupported smoother settings and numerical overflow are checked errors. A valid empty input gives fresh empty
results with preserved metadata and all five empty output arrays; it does not
retain an unrelated old output object as the C++ early-return overload can.
Short or flat traces without a HiRes seed also produce well-formed empty output.

The point limit is checked before smoothing. Gaussian work is bounded before
convolution by its coefficient count plus a conservative count of sample visits
within each kernel's support. SG setup is bounded using frame squared times
polynomial columns, plus sample count times frame length for convolution.
The existing Gaussian coefficient cap and SG numerical limits remain in force.
HiRes work and each active median estimator's work limit are capped by
`max_work`; their own point/bin limits still apply. Boundary lookups, extensions,
overlap walks and all inclusive raw-sum visits share a separate checked work
counter. These are per-stage bounds rather than a claim that every CPU operation
is charged to one global counter. Standard-library binary-search comparisons,
input validation and metadata copies are not separately counted.

Storage is linear in input samples and output peaks, plus the existing smoothing
and noise-estimation tables. No mutable picker cache survives a call. Returned
errors leave an input passed to `filter_chromatogram` unchanged.

Crawdad, transition-group peak selection, multi-transition alignment and
cross-chromatogram consensus are separate remaining ports. This module picks
one chromatogram and does not assign molecular identities or confidence values.

## Evidence

The pinned `PeakPickerChromatogram_test.cpp` supplies two 18-point SRM traces and
four legacy/corrected expected peak, boundary and sample-sum cases. Their source
literals, binary-preserving extracted samples, original assertion precision and
SHA-256 hashes are recorded in
[`chromatogram_processing_provenance.json`](../tests/data/chromatogram_processing_provenance.json).
The independent reference tests retain these source values and verify exact
inclusive sums separately from rounded source assertions.

Focused tests cover irregular sampling, both smoothers, forced-width/noise
interactions, shared valley endpoints, original indices, metadata and omitted
arrays, empty input, validation and atomic numerical failures. Private source
helpers also test rightward nearest ties, the end sentinel and the sequential
overlap midpoint rule. No C++ code was built or executed.
