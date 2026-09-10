# Peak integration

`analysis::peak_integrator::PeakIntegrator` provides read-only spectrum and
chromatogram integration, baseline estimation and peak-shape metrics. It uses
native containers and owned results. By default it operates directly on the
supplied samples. Optional native EMG preprocessing estimates a peak model and
applies the same integration, baseline and shape routines to its fitted samples.

The implementation follows OpenMS4-core revision
[`7c029e8cdba6abab503708ecdd56f6ab55e38ce4`](https://github.com/okohlbacher/OpenMS4-core/tree/7c029e8cdba6abab503708ecdd56f6ab55e38ce4).
The algorithm is in the template implementation of
[`PeakIntegrator.h`](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h).

## API and defaults

`PeakIntegrator` has public typed `integration_method` and `baseline_type`
options, defaulting to `IntensitySum` and `BaseToBase`. The optional
`emg: Option<EmgGradientDescent>` defaults to `None`. The positive `max_points`
limit defaults to 1,000,000 whole-input points and is checked before input
validation or hull allocation; it also caps the number of fitted output points.
Each operation validates the native container,
including finite values, sorted positions and parallel-array lengths. Negative
finite intensities are accepted.

| Operation | Spectrum / chromatogram methods | Result |
| --- | --- | --- |
| Integrate | `integrate_spectrum` / `integrate_chromatogram` | `PeakArea`: area, height, apex, selected `[position, intensity]` hull points |
| Baseline | `estimate_background_spectrum` / `estimate_background_chromatogram` | `PeakBackground`: area and height at the supplied apex |
| Shape | `calculate_shape_metrics_spectrum` / `calculate_shape_metrics_chromatogram` | `PeakShapeMetrics`: widths and sampled starts/ends at 5%, 10%, 50%; total width; tailing/asymmetry; baseline difference/height ratio; sample counts |

Positions are m/z for spectra and retention time in seconds for chromatograms.
Every operation takes finite ordered bounds. Background and shape operations
also require a finite apex inside those requested bounds. Shape height must be
finite and nonnegative. Inputs and their annotations are preserved on success
and error. The API returns checked `Result` values instead of source exceptions,
invalid iterators, or nonfinite calculations.

## Optional EMG preprocessing

Set `emg` to `Some(EmgGradientDescent::default())`, or supply its bounded fitter
configuration. Integration, background and shape operations each fit once using
the requested inclusive bounds, then replace those bounds with the fitted
container's first and last sample. Optional extrapolated tail points therefore
contribute to the entire result, including hull, sampled baseline endpoints,
shape widths and point counts. The caller's supplied background apex and shape
height/apex remain unchanged, as in source `EMGPreProcess_`; they are not replaced
with the fitted maximum. A completely empty input retains the source's special
zero-shape result for valid arguments without fitting. Empty integration or
background fitting requests return checked errors.

The integration bounds are always passed as `Some(left)` and `Some(right)` to
the fitter. A numeric zero is consequently a literal bound. This deliberately
avoids the C++ fitter's zero-means-unbounded sentinel; standalone native fitting
uses `None` to request an unbounded side. Changing the interval fits that subset
of the observed trace, rather than fitting the entire input and cropping later.
Disabling `compute_additional_points` retains the selected sample positions but
still replaces their intensities with fitted values.

The native fitter returns a typed `EmgEstimate` alongside its fitted spectrum or
chromatogram. It does not append the source's four-entry `emg_parameters` array,
which would violate the native peak-alignment invariant. PeakIntegrator's result
shapes stay unchanged; use `EmgGradientDescent::fit_spectrum` or
`fit_chromatogram` directly when parameter estimates, fit diagnostics or omitted
annotation names are needed. Input samples, metadata and annotation arrays are
never mutated. No synthetic annotations are inferred for extrapolated points.
Each operation temporarily owns its fitted container and borrows the original
when fitting is disabled.

The effective fitter point cap is the minimum of `PeakIntegrator::max_points`
and the fitter's own `max_points`, enforced during fitting before extrapolation
can exceed it. The fitter's iteration and evaluation limits also apply. Its
validation and finite-calculation errors propagate unchanged; no failed fit is
silently replaced with an unfitted integration. Signed intensities remain
accepted by the fitter; fitting additionally requires distinct sample positions,
a positive source-initialized mean/width and finite numerical expressions.

## Integration conventions

Bounds select every observed position in the inclusive interval. They do not
interpolate endpoints. Areas do not automatically subtract a baseline.

- `IntensitySum` adds each selected `f32` intensity to an `f64` sum. Its units are
  intensity; the other methods weight intensities by coordinate spacing.
- `Trapezoid` follows the source expression exactly: add neighboring intensities
  in `f32`, promote the sum to `f64`, divide by two, then multiply by position
  difference. Using an `f64` intensity addition changes observable results.
- `Simpson` uses the source three-point formula for nonuniform spacing. With two
  selected points it explicitly falls back to trapezoid. With an odd count of
  at least three, it integrates overlapping three-point segments. With an even
  count of at least four, it averages odd-count subintegrals: omit the last,
  omit the first, include the preceding input sample if available, and include
  the following sample if available, in that order. The latter samples can lie
  outside the requested interval; they do not enter the reported hull or apex.

The even-count source algorithm treats `-1.0` as a missing-subintegral sentinel
and excludes *computed* subareas exactly equal to `-1.0` as well. This unusual
finite behavior is preserved. If all available subareas equal `-1.0`, the native
API returns an error instead of the source's undefined `0/0`. Odd-count Simpson
returns a valid area of `-1.0` normally. Finite negative areas are never clamped;
nonuniform Simpson can yield a negative area even with positive input samples.

Empty selections yield area/height zero and an empty hull. A singleton yields its
intensity under summation and zero under trapezoid/Simpson. The height starts at
zero and updates only for a strictly larger intensity: the first positive
maximum wins ties. Without a positive sample, height stays zero and apex is the
source `(left + right) / 2`. A nonfinite final midpoint returns an error.

## Baseline and shape conventions

`BaselineType` accepts the source names through `FromStr`: `base_to_base`,
`vertical_division_min`, `vertical_division_max`. Legacy `vertical_division`
resolves to `VerticalDivisionMin` and displays using that canonical name.

Base-to-base joins the first/last **sampled** intensities. For summed integration,
the baseline area is the sum of this line at each selected coordinate using the
source rectangle-plus-triangle expression; for trapezoid and Simpson it is the
endpoint trapezoid area. Min/max methods use the lower/higher sampled endpoint
as a constant baseline. Their summed area is baseline times sample count; their
weighted area is baseline times sampled width. Base-to-base endpoint subtraction
occurs in `f64` after both intensities have been promoted. The source height
formula uses the absolute distance from the smaller endpoint; if the supplied
apex is inside the requested bounds but outside the sampled endpoints, this
retains that absolute-distance behavior rather than ordinary line extrapolation.

Shape calculations retain the source sampled threshold search. Starting at each
outer boundary, advance inward while intensity is at or below the requested
fraction, returning the final such sample. If the boundary already exceeds the
threshold, return the boundary. There is no interpolation or global nearest
threshold search; the method assumes a convex peak. The apex partitions the
sample range using the first position at or above it. The right search includes
that sample. Half-height sample counts include equality.

`total_width` is the difference between sampled endpoints. Despite its source
name, `slope_of_baseline` is the right-minus-left **intensity difference**, without
dividing by width; that subtraction occurs in `f32` before promotion. Tailing is
width-at-5% divided by twice the apex-to-start distance. Asymmetry uses the
end-to-apex divided by apex-to-start distances at 10%. Source guarded conventions
are retained: tailing/asymmetry are zero when the respective start equals the
apex; the baseline-to-height ratio is zero when supplied height is zero. Those
zeros can indicate a degenerate sampled peak and do not imply an ideal shape.

## Checked boundaries and remaining scope

An empty whole trace returns zero shape metrics for valid arguments. A nonempty
trace with no selected sample, or a shape apex whose lower bound lies beyond the
selected samples, returns an error. This avoids reading a sample outside the
requested interval. A singleton at its apex has zero widths and guarded shape
ratios. Empty background selections are errors. Base-to-base requires distinct
sampled endpoint positions; min/max baselines remain defined for singletons or
coincident endpoints (zero weighted area, ordinary summed area).

Repeated positions are accepted for summation and trapezoid, and for shape
metrics when their guarded expressions are defined. Any repeated spacing used
in a Simpson subintegral is an error. Nonfinite results and overflowing tailing
factor denominators are errors. The implementation does not silently change the
formula to compensate for overflow or underflow.

There is no iterator-overload framework, automatic background subtraction, or
background-corrected shape calculation. Without fitting, work in the scientific
routines is linear in whole-input point count; integration hull storage is linear
in selected points, and baseline/shape routines allocate no point arrays. With
fitting, additional work and storage follow the fitter limits and include the
owned fitted container. Native container validation also checks any attached
metadata/identification records, which have their own size and validation costs.

## Validation and provenance

`tests/peak_integrator.rs` checks typed defaults and aliases, inclusive bounds,
input preservation, signed and empty/singleton behavior, baseline choices and
directions, clipped threshold searches, duplicates, invalid arguments, resource
limits and arithmetic overflow. Independent
`tests/chromatogram_processing_reference.rs` uses the exact 107-point upstream
L-glutamate trace and the three-point negative Simpson regression. It verifies
source integration/baseline/shape goldens, even-count neighboring samples,
float-intensity precision and both sentinel collision outcomes. The fixture
provenance records literal source assertions and independent analytical cases;
no C++ program was built or executed. `tests/emg_integration.rs` additionally
checks all three integration methods and all baseline choices after fitting,
expanded-span shape calculations with unchanged supplied height/apex, cropped
inputs, zero-boundary semantics, no-fit defaults, bounded extrapolation, and
input/annotation preservation on both success and error. Its cutoff fixtures
restore 12→28 points on the right and 66→71 on the left, using exact source bits
from `tests/data/emg_traces.tsv`. The expected fitted lengths come directly from
the pinned `EmgGradientDescent_test.cpp`; fitter-specific parameter and curve
oracles are tested separately.

Pinned SHA-256 values:

| Source file | SHA-256 |
| --- | --- |
| `src/openms/include/OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h` | `93fbfa46f512181bd8a646026014f56a9614222dd1222cecb6fc8785ea693c86` |
| `src/openms/source/ANALYSIS/OPENSWATH/PeakIntegrator.cpp` | `0c88cb7d0c88dd618039ae641186993af958f8d4b41e9bbbd8662bcba70367fb` |
| `src/tests/class_tests/openms/source/PeakIntegrator_test.cpp` | `dfda657e349d84ec132c60a3ab8e554178e9a485d3d4034758b1eed3add56a2f` |
