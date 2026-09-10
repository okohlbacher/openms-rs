# Retention-time transformation models

`analysis::transformations` provides native Rust implementations of `TransformationDescription`, linear regression with coordinate weighting, linear and natural cubic interpolation, robust LOWESS followed by interpolation, and application to the represented map and identification containers. The reference is OpenMS4-core commit `7c029e8cdba6abab503708ecdd56f6ab55e38ce4`. The implementation calls no C++ library and adds no dependencies. This is a completed set of model implementations within the ongoing core port; it does not implement the complete map-alignment subsystem.

```rust
use openms::analysis::transformations::{
    DataPoint, ModelConfig, TransformationDescription,
};

fn aligned_retention_times(values: &mut [f64]) -> openms::Result<()> {
    let mut transformation = TransformationDescription::new(vec![
        DataPoint::with_note(10.0, 12.0, "landmark A"),
        DataPoint::with_note(20.0, 23.0, "landmark B"),
        DataPoint::with_note(30.0, 31.0, "landmark C"),
    ])?;
    transformation.fit_model(ModelConfig::Lowess(Default::default()))?;
    transformation.apply_values(values)
}
```

Coordinates are `f64` values. The models have no implicit time-unit conversion: anchors and queries must use the intended source and target units. Negative finite coordinates and responses are valid unless a chosen coordinate weight makes them undefined.

## Models and source behavior

| Native model | Implemented behavior |
| --- | --- |
| `None` | Passes coordinates through; this is the initial unfitted state. |
| `Identity` | Passes coordinates through and ignores subsequent fit requests until anchors are replaced, following the source identity lock. |
| `LinearModel` | Ordinary least squares with an intercept; one anchor gives slope one and the observed offset, two anchors use their exact line, and three or more use centered covariance and variance. Explicit slope/intercept are accepted only as the model when no anchors are supplied. |
| `InterpolatedModel` | Sorts anchors, averages responses at repeated source coordinates, and uses linear or natural cubic interpolation. At least three distinct source coordinates are required, including for the linear option. |
| `LowessModel` | Runs the source robust locally weighted linear smoother and then the selected interpolation model on its fitted observations. The public `lowess` helper also returns smoothed observations directly. |

`CoordinateWeight` exposes identity, natural logarithm, reciprocal absolute value, and reciprocal squared value. These transform the x and y coordinates before regression; they are not observation-specific regression weights. Training values are clamped to the configured range only when that coordinate's weight is nonidentity. Default bounds are `1e-15` and `1e15`. Evaluation deliberately does not clamp coordinates. The response transform is undone after evaluating the fitted line. Reciprocal transforms discard signs, and their inverses describe the nonnegative branch.

Interpolation uses the existing native `CubicSpline2d` natural boundary conditions and tridiagonal recurrence. It does not use a smoothing spline. All three source extrapolation policies are available:

- `TwoPointLinear`: one line through the first and last averaged anchor; the interpolation model's default.
- `FourPointLinear`: separate lines through the first two and last two averaged anchors; the LOWESS model's default.
- `GlobalLinear`: ordinary least squares over every original anchor, retaining the contribution of repeated source coordinates. In LOWESS, these are the smoothed observations before duplicate averaging.

LOWESS defaults are span `2/3`, three robustifying iterations after the initial fit, and delta equal to one percent of the coordinate range. `Some(delta)` supplies an explicit finite nonnegative distance; `None` selects automatic delta. The neighborhood contains `floor(span*n)` observations, clamped to the interval from two through n. The port retains tricube distance weights, local linear correction, robustness reweighting, ties, and linear interpolation over delta-skipped fits. It also retains the pinned implementation's residual pseudo-median calculation, including its two selected central residuals for an odd observation count. This detail is needed to match the actual implementation rather than a different LOWESS variant.

The public `lowess` helper requires at least two observations sorted by nondecreasing x. `LowessModel::fit` sorts its anchors and subsequently needs three distinct x values for interpolation. Signed responses are supported. Rust preserves input order among equal coordinates; C++ uses an unstable sort with a comparator on x alone. Minor rounding differences or robust tie-order effects can therefore occur for duplicate coordinates even though the source fixture cases agree.

## Descriptions, inverses, and diagnostics

`TransformationDescription` keeps the original annotated anchors separate from its fitted model. Replacing anchors resets the model to `None`; failed replacements or fits leave the previous state intact. `apply` evaluates one coordinate, and `apply_values` applies a batch atomically. Notes stay attached to their anchor when it is inverted.

`TransformationDescription::inverse` swaps the source and target anchors and refits the same model configuration, as in OpenMS. `invert` replaces the description only after this succeeds. This is usually an approximate inverse: noisy forward and backward least-squares fits are different, and an interpolated curve can be nonmonotonic. The implementation does not solve the original curve's inverse or guarantee monotonicity. Inversion can fail when swapped coordinates do not support the selected model. Coordinate weights stay in their original configuration when an anchored description is refitted, matching the source.

`LinearModel::inverse`, and description inversion for a linear model with no anchors, use the algebraic reciprocal slope and adjusted intercept and exchange the x/y weight configurations. A zero slope is an error.

`deviations` reports absolute anchor errors before or after transformation, optionally sorted. `statistics` returns source/target ranges and the source percentile selections at 100, 99, 95, 90, 75, 50, and 25 percent. These statistics use the source discrete rank formula rather than interpolated quantiles. For small samples, Rust clamps the resulting negative rank to the first observation. Empty input returns absent ranges and no percentile records.

`estimate_window` follows the source adaptive residual quantile: compute the requested interpolated quantile, cap observations at `Q3 + 1.5*IQR` for the robust estimate, and blend robust and uncapped estimates as tail density rises from one to ten percent. Fewer than four residuals or zero IQR use the uncapped estimate. Defaults are quantile 0.99, inversion before residual calculation, full window width rather than half width, and padding factor one. Empty residuals produce zero. The source-unit/default window uses the refitted inverse described above.

## Applying transformations to data

`analysis::alignment_transformer::MapAlignmentTransformer` applies a description to experiments, feature maps, consensus maps, or slices of peptide identifications. It also exposes single-feature and single-consensus-feature methods. Every operation prepares its changes before mutating the input; a failure in the final chromatogram sample, subordinate hull, consensus handle, or unassigned identification leaves the entire supplied container unchanged.

```rust
use openms::analysis::alignment_transformer::MapAlignmentTransformer;
use openms::analysis::transformations::TransformationDescription;
use openms::MSExperiment;

fn align(experiment: &mut MSExperiment, transformation: &TransformationDescription)
    -> openms::Result<()>
{
    MapAlignmentTransformer {
        store_original_rt: true,
        ..Default::default()
    }.transform_experiment(experiment, transformation)
}
```

| Container | Retention times transformed |
| --- | --- |
| Experiment | Every spectrum RT and every chromatogram sample RT. |
| Feature map | Feature centroids, each stored hull outline point, all subordinate features, assigned peptide IDs, and unassigned peptide IDs. |
| Consensus map | Consensus centroids, every source-feature handle, assigned peptide IDs, and unassigned peptide IDs. Centroids are transformed directly, without recomputing a consensus from its handles. |
| Peptide identification slice | Existing `Some(rt)` values; records without an RT remain unchanged. |

The source experiment overload leaves attached spectrum peptide identification RTs unchanged. The Rust default preserves that behavior; `transform_spectrum_identifications: true` enables their transformation as an explicit extension. Feature and consensus identification RTs always transform. Protein identification records, evidence, scores, masses, intensities, array values, identifiers, and other metadata retain their values.

When `store_original_rt` is enabled, scalar records gain `original_RT` only if that key is absent. Existing values are retained on repeated transformations even when their type or contents differ from a newly generated value. Peptide, feature and consensus metadata stores typed floating values. Spectrum metadata still uses string maps, so spectrum scalar values use round-trippable decimal strings. Chromatograms use the source's lower-case key `original_rt` and preserve the entire original time vector as a string such as `[1, 3]`; empty chromatograms store `[]`. This string encoding reflects the current kernel metadata representation and is not a typed C++ `DataValue` list.

The transformer keeps spectrum, feature, chromatogram-sample, and annotation order. Nonmonotonic or decreasing transformations can make retention-time sequences unsorted; call the kernel's sorting operations explicitly when a downstream algorithm requires sorted times. Aligned data arrays keep their original sample association. Native ranges are computed from current data and do not use stale cached bounds.

For hulls, the source reads its ordered outline, transforms the RT coordinate of each vertex, and stores the transformed outline without rebuilding the geometric convex hull. The port follows this behavior even for nonmonotonic transformations. It therefore converts a scan-envelope hull to an outline-only hull. Current kernel scan-envelope containment queries reject outline-only hulls; bounding boxes and outline vertices remain available. The transformer does not infer an inverse, restore lost scan envelopes, resample hull edges, or correct self-intersections introduced by a model.

Transformer defaults limit each operation to one million transformed coordinates, one million spectra/chromatograms/features/hulls/peptide records, and subordinate depth 128. Hull vertices and consensus handles count as transformed coordinates. Missing-RT peptide records count toward the record limit. Limits are configurable, with depth capped by the kernel limit. Input validation also traverses peak arrays and identification contents; the limits do not cap their byte sizes, annotation strings, or the cost of that validation. Edit plans store transformed RTs and hulls rather than copying spectral peak data. Consensus features are copied during preparation. Existing hull outlines must be materialized before their coordinate count is checked. Allocation failure follows the allocator's normal behavior.

The broader `IdentificationData` observation model and its transformer overload are outside the current port. These methods apply existing transformations; they do not calculate map correspondences.

## Validation and resource limits

Anchor coordinates, queries, coefficients, bounds, and arithmetic results must be finite. Singular linear fits, log/reciprocal domain failures, unrepresentable coordinate differences, insufficient distinct interpolation coordinates, and arithmetic overflow return errors. LOWESS span must be greater than zero and at most one, delta must be finite and nonnegative when explicit, and robustifying iterations are capped at 64. Invalid parameter states cannot silently select another algorithm.

Each model defaults to a maximum of one million observations; the description separately limits stored anchors and batch size. LOWESS additionally defaults to fifty million counted local sample, interpolation, and residual operations. This is a deterministic loop-work guard, not a wall-clock limit. Sorting is bounded by observation and iteration limits but is not included in that work counter. Models and fitting use linear auxiliary storage, with multiple owned arrays; they do not stream large anchor sets. Allocation failure remains subject to the Rust allocator's normal process behavior.

Finite checks and explicit limits intentionally improve undefined or nonfinite source failure cases. For example, nonfinite residuals produce an error instead of being discarded from a window estimate. Fit, anchor replacement, in-place inversion, and batch application do not partially mutate their inputs on a reported error. Floating-point underflow and ordinary rounding are not treated as errors.

## Remaining scope

Akima interpolation, B-spline fitting, LOWESS automatic span selection/cross-validation and its scoring/grid helpers, TrafoXML I/O, and the C++ string-based `Param` API are not implemented here. No natural-spline approximation is presented as the source B-spline model. The C++ `symmetric_regression` flag is assigned but unused in the pinned linear implementation; the Rust API exposes the actual ordinary fit and has no ineffective flag. Native interpolation follows the primary annotated-`DataPoint` overload and does not reproduce alternate C++ container-overload allocation defects.

These models fit supplied anchor correspondences. They do not discover landmarks, match maps, or implement the broader alignment algorithms. See the port status document for coverage of the surrounding kernel and processing modules.

## Verification and provenance

`tests/transformations.rs` checks the pinned weighted/unweighted linear examples, exact and refitted inverses, duplicate-anchor averaging, all extrapolation policies, the upstream interpolation table, LOWESS sine predictions, both fifty-observation cars fits, and the original twenty-observation LOWESS delta/robustness examples. It also checks notes, identity locking, source deviations and percentiles, adaptive windows, failed-operation atomicity, finite-domain rejection, and resource limits. `tests/alignment_transformer.rs` reuses the four container examples from the source transformer test and checks preservation of original RTs, chromatogram arrays, decreasing transforms, subordinate hulls, consensus handles, assigned/unassigned peptide RTs, the optional spectrum-ID extension, limits, and late-error atomicity.

The fixture manifest [transformations_provenance.json](../tests/data/transformations_provenance.json) records source and fixture SHA-256 hashes and extraction details. Interpolation tables are byte-for-byte copies; LOWESS fixtures extract numeric literals from upstream tests. No expected values were generated by the Rust implementation and no C++ build was performed. Tolerances reflect the precision of each source fixture: the interpolated table uses absolute `1e-5`, the cars literals `1e-10`, the rounded original LOWESS examples `1e-3`, and sine predictions `1e-6`.

The port retains the OpenMS BSD-3-Clause header. The pinned `FastLowessSmoothing` source credits W. S. Cleveland's NETLIB LOWESS and a BSD translation in common-lisp-stat; this provenance is retained in the Rust smoother. No third-party dependency source is modified.

Pinned sources:

- [TransformationDescription](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationDescription.cpp) and [tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/TransformationDescription_test.cpp).
- [Linear model](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelLinear.cpp), [coordinate weights](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModel.cpp), and [linear tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/TransformationModelLinear_test.cpp).
- [Interpolated model](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.cpp) and [tests](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/TransformationModelInterpolated_test.cpp).
- [LOWESS model](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/TransformationModelLowess.cpp), [FastLowessSmoothing](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/PROCESSING/SMOOTHING/FastLowessSmoothing.cpp), and [numerical examples](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/FastLowessSmoothing_test.cpp).
- [Adaptive quantile and quantile helpers](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/include/OpenMS/MATH/StatisticFunctions.h).

- [MapAlignmentTransformer](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/openms/source/ANALYSIS/MAPMATCHING/MapAlignmentTransformer.cpp) and [container examples](https://github.com/okohlbacher/OpenMS4-core/blob/7c029e8cdba6abab503708ecdd56f6ab55e38ce4/src/tests/class_tests/openms/source/MapAlignmentTransformer_test.cpp).
