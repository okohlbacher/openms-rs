# Retention-time transformation models

`analysis::transformations` ports the retention-time transformation headers of `ANALYSIS/MAPMATCHING`: `TransformationDescription`, `TransformationModel`, `TransformationModelLinear`, `TransformationModelInterpolated`, `TransformationModelBSpline`, `TransformationModelLowess` and `TransformationModelDefaults`, plus the application of a transformation to maps (`MapAlignmentTransformer`, in `analysis::alignment_transformer`). The reference is OpenMS4-core `bc9cc12514c768385ce121d6ca4bb710fe1983c4`; the earliest fixtures were extracted at `7c029e8`, and the ones that changed since are re-pinned in [transformation_models_provenance.json](../tests/data/transformation_models_provenance.json). The implementation calls no C++ library and adds no dependency.

| Header | Rust | Status |
|---|---|---|
| `TransformationDescription.h` | `src/analysis/transformations.rs` | complete (this document) |
| `TransformationModel.h` | `src/analysis/transformations.rs`, `transformations/params.rs` | complete |
| `TransformationModelLinear.h` | `src/analysis/transformations.rs` | complete |
| `TransformationModelInterpolated.h` | `src/analysis/transformations.rs`, `transformations/akima.rs` | complete |
| `TransformationModelBSpline.h` | `src/analysis/transformations/b_spline.rs` | complete |
| `TransformationModelLowess.h` | `src/analysis/transformations/lowess.rs` | complete |
| `TransformationModelDefaults.h` | `src/analysis/transformations/params.rs` | complete |
| `MapAlignmentTransformer.h` | `src/analysis/alignment_transformer.rs` | complete |

```rust
use openms::analysis::transformations::{DataPoint, ModelConfig, TransformationDescription};
use openms::param::Param;

fn aligned_retention_times(values: &mut [f64]) -> openms::Result<()> {
    let mut transformation = TransformationDescription::new(vec![
        DataPoint::with_note(10.0, 12.0, "landmark A"),
        DataPoint::with_note(20.0, 23.0, "landmark B"),
        DataPoint::with_note(30.0, 31.0, "landmark C"),
    ])?;
    // Typed, or by name from a Param exactly as `fitModel(model_type, params)`:
    transformation.fit_model(ModelConfig::Lowess(Default::default()))?;
    transformation.fit_model_with_param("linear", &Param::new())?;
    transformation.apply_values(values)
}
```

## Two interfaces to one model family

The source's models are a virtual base class with a subclass per model, each constructed from `(DataPoints, Param)`. The port keeps the models as a closed enum, `TransformationModel`, fitted from a typed `ModelConfig` whose option structs (`LinearOptions`, `InterpolationOptions`, `BSplineOptions`, `LowessOptions` with `LowessAutoSpan`) are the primary API. The source's `Param` interface is kept beside it, because the alignment tools pass model parameters through by name:

- every options struct has `defaults() -> Result<Param>` (the class's static `getDefaultParameters`, same names, values, types, descriptions and restrictions, in declaration order), `from_param(&Param)` (the constructor's reading of its `Param`) and `to_param()`;
- `ModelConfig::from_param(model_type, &Param)` is the model choice `fitModel` makes, and `ModelConfig::default_parameters(model_type)` the matching defaults;
- `TransformationDescription::fit_model_with_param(model_type, &Param)` is `fitModel(model_type, params)`, and `model_parameters()` returns exactly what the source model keeps in `params_` (see below).

A model constructor in the source is not a `DefaultParamHandler`: it fills missing values from its defaults (`Param::setDefaults`), checks no restriction and warns about no unknown name. `from_param` does the same. It refuses only what the source's `ParamValue` conversion refuses (`Exception::ConversionError`: a string where a number is read, a negative count for an unsigned one) and what the constructor rejects itself.

## API mapping

### `TransformationDescription`

| Source | Rust |
|---|---|
| `DataPoint`, `DataPoints` (typedefs of the model's) | `DataPoint`, `Vec<DataPoint>` |
| `TransformationStatistics` (`percents`, `xmin`..`ymax`, `percentiles_before/after`) | `TransformationStatistics` (`x_range`, `y_range` as `Option<NumericRange>`, `percentiles: Vec<DeviationPercentile>` in the source's 100, 99, 95, 90, 75, 50, 25 order) |
| default constructor, constructor from `DataPoints` | `Default`, `new`, `with_max_points` (native limit) |
| copy constructor, `operator=` | `Clone` (copies the fitted model; see native differences) |
| `fitModel(model_type, params = Param())` | `fit_model_with_param(model_type, &params)`; typed `fit_model(ModelConfig)` |
| `apply` | `apply`; batch `apply_values` (native, atomic) |
| `getModelType` | `model_type()` (and `model().name()`) |
| static `getModelTypes(StringList&)` | `model_types()` / `MODEL_TYPES`: `linear`, `b_spline`, `interpolated`, `lowess` |
| `setDataPoints(DataPoints)` | `set_data_points` |
| `setDataPoints(vector<pair<double,double>>)` | `set_data_point_pairs` |
| `getDataPoints` | `data_points` |
| `getModelParameters` | `model_parameters` |
| `invert` | `invert` (atomic), non-mutating `inverse` |
| `getDeviations(diffs, do_apply, do_sort)` | `deviations(apply, sort)` |
| `getStatistics` | `statistics` |
| `printSummary(ostream&)` | `print_summary(&mut impl Write)` |
| `estimateWindow(quantile, invert, full_window, padding_factor)` | `estimate_window(WindowOptions)` with the same defaults |
| protected `data_`, `model_type_`, `model_` | private fields; `model()` and `model_config()` read them |

### `TransformationModel`

| Source | Rust |
|---|---|
| `DataPoint` (`first`, `second`, `note`; pair constructor; `operator<` and `operator==` over `std::tie`) | `DataPoint` (`x`, `y`, `note`; `From<(f64, f64)>`; derived `PartialOrd`/`PartialEq` in the same field order) |
| default constructor (identity) | `TransformationModel::None` / `Identity` |
| constructor `(DataPoints, Param)` | the base reads only weights: `ModelWeights::from_param` |
| `evaluate` | `TransformationModel::apply` and each model's `apply` |
| `weightData`, `unWeightData` | `ModelWeights::weight_data` / `unweight_data`; `LinearModel::weight_data`; `TransformationModel::weight_data` (a no-op for the unweighted models, see below) |
| `checkValidWeight` | `check_valid_weight` |
| `checkDatumRange` | `check_datum_range` |
| `weightDatum`, `unWeightDatum` (by name) | `weight_datum`, `unweight_datum`; typed `CoordinateWeight::transform` / `untransform` |
| `getParameters` | `TransformationDescription::model_parameters`; `ModelConfig::to_param`. A base model built from a `Param` reports that `Param` unchanged; `ModelWeights` does not keep a copy, the caller holds it |
| `getValidXWeights`, `getValidYWeights` | `VALID_X_WEIGHTS`, `VALID_Y_WEIGHTS` |
| static `getDefaultParameters` (clears) | `ModelConfig::default_parameters("none")` (empty) |
| protected `x_weight_`, `x_datum_min_`.. `weighting_` | `ModelWeights { x, y }` of `CoordinateWeight`; `is_weighting()` |

### `TransformationModelLinear`

Complete since earlier; see "Models and source behaviour". `getDefaultParameters` is `LinearOptions::defaults()` (including the unused `symmetric_regression`), and `getParameters(slope, intercept, ...)` is `coefficients()` plus `options()`.

### `TransformationModelInterpolated`

| Source | Rust |
|---|---|
| constructor `(DataPoints, Param)` | `InterpolatedModel::fit(data, InterpolationOptions)`; `InterpolationOptions::from_param` |
| constructor `(vector<pair>, Param, bool preprocess)` | `InterpolatedModel::fit_pairs(pairs, options, preprocess)`, reproducing its zero padding (see below) |
| `evaluate` | `apply` |
| static `getDefaultParameters` | `InterpolationOptions::defaults()` |
| nested abstract `Interpolator` (`init`, `eval`) and the `.cpp`'s `LinearInterpolator`, `Spline2dInterpolator`, `AkimaInterpolator` | the closed `Interpolation` enum: `Linear`, `CubicSpline` (the crate's `CubicSpline2d`), `Akima` (`akima.rs`, a port of GTE `IntpAkimaNonuniform1`) |
| private `preprocessDataPoints_` (both overloads) | private `averaged` |

### `TransformationModelBSpline`

| Source | Rust |
|---|---|
| constructor `(DataPoints, Param)` | `BSplineModel::fit(data, BSplineOptions)`; `BSplineOptions::from_param` |
| `evaluate` | `apply` |
| `getParameters` | `model_parameters` of the description |
| static `getDefaultParameters` | `BSplineOptions::defaults()` |
| protected `spline_`, `xmin_`/`xmax_`, `extrapolate_` (`EX_LINEAR`, `EX_BSPLINE`, `EX_CONSTANT`, `EX_GLOBAL_LINEAR`), offsets and slopes | `spline()`, `data_range()`, `BSplineExtrapolation`, private fields |

### `TransformationModelLowess`

| Source | Rust |
|---|---|
| constructor `(DataPoints, Param)` | `LowessModel::fit` (fixed span) and `LowessModel::fit_auto` (`auto_span`); `LowessOptions::from_param` returns the options and the `LowessAutoSpan` when enabled |
| `evaluate` | `apply` |
| `getParameters` | `model_parameters` of the description |
| static `getDefaultParameters` | `LowessOptions::defaults()` (all twelve parameters, the automatic-selection ones included) |
| private `buildSpanGrid` | `build_span_grid` |
| private `CVMetric`, `names_of_cvmetric` | `CrossValidationMetric`, `CrossValidationMetric::ALL` and `name`/`from_name` |
| private `scoreResiduals` | `score_residuals` |
| private `kTieTol` | `SPAN_SCORE_TIE_TOLERANCE` |
| the grid search (via `CrossValidation::makeKFolds` and `gridSearch1D`) | `LowessModel::select_span` |
| protected `model_` | the interpolation inside `LowessModel` (`knots()`) |
| `FastLowessSmoothing::lowess` | `lowess` |

### `TransformationModelDefaults`

| Source | Rust |
|---|---|
| static `getDefaults(default_model)` | `transformation_model_defaults(default_model)` |

## Models and source behaviour

| Native model | Implemented behaviour |
| --- | --- |
| `None` | Passes coordinates through; the initial, unfitted state. |
| `Identity` | Passes coordinates through and ignores later fit requests until the anchors are replaced, the source's identity lock. |
| `LinearModel` | Ordinary least squares with an intercept; one anchor gives slope one and the observed offset, two anchors their exact line, three or more centred covariance over variance. Explicit slope and intercept are the model only when no anchors are supplied. |
| `InterpolatedModel` | Sorts anchors, averages responses at repeated x (in anchor order, from 0.0 as `std::accumulate` does), and interpolates linearly (by `std::upper_bound`), with a natural cubic spline, or with an Akima spline. At least three distinct x are required, also for linear. |
| `BSplineModel` | The crate's `BSpline2d` (Ooyama's smoothing spline, `docs/BSPLINE2D_SUPPORT.md`) with the model's checks and four extrapolations: `linear` (the spline's slope at the end), `b_spline`, `constant` (the spline's value at the end), `global_linear` (a least-squares line through all anchors, evaluated at the ends). |
| `LowessModel` | Sorts the anchors, runs the robust local linear smoother, and interpolates the smoothed points with the chosen interpolation and extrapolation; optionally selects the span by cross-validation first. |

`CoordinateWeight` exposes identity, natural logarithm, reciprocal absolute value and reciprocal square. These transform x and y before regression; they are not per-observation regression weights. Training values are clamped to the datum range only for a nonidentity weight; the default bounds are `1e-15` and `1e15`; evaluation never clamps. The response transform is undone after evaluating the line. Reciprocal transforms discard signs, and their inverses give the non-negative branch.

**B-spline.** Fewer than two anchors, a `wavelength` larger than the x range, and a spline that cannot be fitted are errors, as the source's `IllegalArgument`, `InvalidValue` and `UnableToFit`. An `extrapolate` value other than `b_spline`, `global_linear` and `constant` means linear, as in the source, which compares only those three; only the parameter's valid-string restriction names the four values.

**Akima.** The derivative at each knot is Akima's weighted average of the neighbouring chord slopes, with two slopes extrapolated linearly at each end (`2 s2 - s3`, ...), and each interval is the Hermite cubic through the knot values and derivatives. Evaluation clamps to the knot range first and finds the interval by a linear scan, as GTE does; the model only evaluates inside the range, so the clamp is never reached there.

**LOWESS automatic span.** With `auto_span`, the anchors are sorted, the candidate grid built by `build_span_grid` (the default grid 0.01 .. 0.90, clamped to `[max(0.01, auto_span_min), min(0.99, auto_span_max)]` with the lower end raised to `min(0.9, auto_min_neighbors / n)`, sorted, near-duplicates within `1e-9` removed), and every candidate cross-validated: leave-one-out for up to 50 anchors, `max(2, auto_k_folds)` round-robin folds otherwise; a fold is skipped when fewer than `max(3, auto_min_neighbors)` anchors remain; a candidate with `span * n < auto_min_neighbors` scores `+inf`. Scores are `rmse`, `mae` (the default) or the type-7 90/95/99th percentile of the absolute held-out errors; a later candidate wins by more than `1e-12`, or by a tie and a larger span. The selected span is fitted and persisted: the stored parameters carry it with `auto_span` set to `false`, and the description's `model_config()` is a fixed-span `Lowess`, so a refit or an inversion does not select again.

**Parameters the model keeps.** `model_parameters()` after `fit_model_with_param` is what the source model's `params_` holds: nothing for `none` and `identity` (the base default constructor ignores its `Param`); for `linear` with explicit `slope` and `intercept` and no anchors, the given parameters unchanged; for `linear` otherwise, the given parameters with the linear defaults filled in and the fitted `slope` and `intercept` set; for the other models, the given parameters with the model's defaults filled in, plus the selected span for an automatic LOWESS. After the typed `fit_model`, the configuration's `to_param` in the same layout. `invert` updates them as the source does: a refitted linear model sets the new `slope` and `intercept`, and the algebraic inverse of an explicit linear model also swaps the weight names and all four datum bounds.

**Statistics and summary.** `statistics` returns the x and y ranges and the percentile selections at 100, 99, 95, 90, 75, 50 and 25 percent, taken at rank `p / 100 * n - 1` of the sorted deviations, the source's discrete formula. `print_summary` writes the source's text exactly, numbers at six significant digits as `std::ostream` prints them, and the "after" block only when a model other than `none` or `identity` is fitted.

`estimate_window` follows the source's adaptive residual quantile: the requested interpolated quantile, the robust quantile of the residuals capped at `Q3 + 1.5 * IQR`, and a blend of the two as the tail beyond the cap grows from one to ten percent. Fewer than four residuals or a zero IQR use the uncapped quantile. The defaults are quantile 0.99, residuals in original units through the refitted inverse, full width and padding one; no anchors give zero.

## Native differences

- **Cloning copies the model.** The source's copy constructor and assignment refit the model from `getModelParameters()`. That gives the same model where it succeeds, but it throws for an inverted explicit linear model: the inversion stores the old y weight name (`y`) as the x weight, which the linear constructor's weight check then rejects. A clone here cannot fail.
- **Atomic inversion.** The source's `invert` swaps the anchors before refitting, so a failed refit leaves them swapped; `invert` here changes nothing on error.
- **Non-finite values are refused.** Coordinates, queries, coefficients, bounds and results must be finite; `ln(0)` in a weighted fit, a zero-variance regression and a spline over a zero interval return errors where the source computes `-inf` or NaN. The string helpers `weight_datum` and `unweight_datum` are the exception: they reproduce the source's arithmetic, `-inf` included, because the class test asserts it.
- **`weightData` on unweighted models.** Only the linear model and the base `(DataPoints, Param)` constructor initialise the weight members; the interpolated, B-spline and LOWESS models are built with the base default constructor, so their `weightData` reads an uninitialised flag. `TransformationModel::weight_data` leaves the data unchanged for them, which is what their `evaluate` implies.
- **Sort order of LOWESS anchors.** The anchors are sorted with the permutation the Release build's `std::sort` produces under the source's x-only comparator (`math::source_sort`), so anchors with equal x reach the smoother in the source's order. Before this lane the port sorted stably and documented the difference.
- **Resource limits.** Each model defaults to at most one million anchors; LOWESS counts at most fifty million local operations per fit; the automatic span selection refuses more than `LowessAutoSpan::max_fits` (10,000) cross-validation fits before the first; the B-spline inherits `BSpline2d`'s point and node ceilings. None of these exists in the source.
- **Stricter than the source, by its own preconditions.** A LOWESS `span` of zero with `auto_span` off is refused: `FastLowessSmoothing::lowess` states `0 < f <= 1` as a precondition checked only in debug builds. A negative `num_nodes` or `num_iterations`, and a `boundary_condition` outside 0..2 (which the source casts to its enum unchecked), are refused.
- **The pair constructor's zero padding is reproduced.** `TransformationModelInterpolated(vector<pair>, Param, preprocess)` resizes its vectors and then appends to them: without `preprocess` the interpolation runs on as many `(0, 0)` points as there are pairs followed by the pairs; with it, `global-linear` extrapolation fits its line to that many `(0, 0)` points and the averaged anchors. `fit_pairs` reproduces both, deterministic as they are, and refuses what the source cannot compute: a cubic or Akima spline over the repeated zero abscissae. Recorded as a defect candidate in the lane report.
- **`fitModel` with an unknown `auto_metric`** fails even with `auto_span` off, as in the source, which converts it on every fit.

## Verification and provenance

`tests/transformation_models.rs` (17 tests) checks, from the pinned class tests: the 41 B-spline predictions at `v = -4, -3.8, ...` (the loop variable accumulated as the C++ does) and the linear, constant and global-linear extrapolation pairs; the B-spline's two-point minimum and its reported parameters; `fitModel` by name for linear, b_spline, lowess and identity on the `TransformationDescription` data, including the `[EXTRA]` section where one `Param` carries a LOWESS and a B-spline setting; `getModelParameters` (empty, and the eight explicit linear parameters returned unchanged); the exact `printSummary` text before and after a linear fit; `getModelTypes`; `TransformationModelDefaults` for six selected models; the LOWESS auto-span tie section (span 0.8 persisted, `auto_span` false, `f(0) = 1`, `f(5) = 11`, `f(10) = 21`); the LOWESS sine predictions through `Param`; Akima against the GSL table's third column; the interpolated defaults and extrapolation through `Param`; and the base model's weight helpers and `DataPoint` ordering. Derived checks cover `build_span_grid` and `score_residuals`, the pair constructor's padding (least squares over six points computed by hand), and the inverted explicit linear parameters. Tolerances are ClassTest's `isRealSimilar` (absolute `1e-5` or ratio `1 + 1e-5`).

`tests/transformations.rs` keeps the earlier checks: weighted and unweighted linear examples, exact and refitted inverses, duplicate-anchor averaging, all extrapolation policies, the interpolation table, LOWESS sine predictions, both fifty-observation cars fits, the original twenty-observation LOWESS delta and robustness examples, notes, identity locking, deviations and percentiles, adaptive windows, atomicity, finite-domain rejection and limits. `tests/transformation_xml.rs` now reads a `b_spline` record to the B-spline defaults (it asserted `Unsupported` while the model was missing) and accepts `akima`.

The fixture manifests are [transformations_provenance.json](../tests/data/transformations_provenance.json) (earlier fixtures, pin `7c029e8`) and [transformation_models_provenance.json](../tests/data/transformation_models_provenance.json) (pin `bc9cc12`, including the hashes of the vendored GTE Akima headers the port was read from). No expected value was generated by the Rust implementation and no C++ was executed for these models; the B-spline numbers are independently backed by `BSpline2d`'s executed probe (`docs/BSPLINE2D_SUPPORT.md`), which evaluates the same 41 points.

The port retains the OpenMS BSD-3-Clause header. The pinned `FastLowessSmoothing` source credits W. S. Cleveland's NETLIB LOWESS and a BSD translation in common-lisp-stat, and `akima.rs` credits the Geometric Tools Engine (Boost Software License 1.0); both attributions are kept in the Rust sources.

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

`transform_identification_data` is the `IdentificationData` overload: every observation with an RT is transformed, with `original_RT` stored under the same rule, in one graph transaction; an observation without RT is left alone, where the source's NaN RT would be transformed to NaN and stored as the original. Molecules, matches and scores are untouched. These methods apply existing transformations; they do not calculate map correspondences — the alignment algorithms that do are in `analysis::mapmatching` ([MAP_ALIGNMENT_SUPPORT.md](MAP_ALIGNMENT_SUPPORT.md)).

| `MapAlignmentTransformer` (source) | Rust |
|---|---|
| `transformRetentionTimes(PeakMap&, trafo, store_original_rt)` | `transform_experiment` |
| `transformRetentionTimes(FeatureMap&, ...)` | `transform_feature_map` |
| `transformRetentionTimes(ConsensusMap&, ...)` | `transform_consensus_map` |
| `transformRetentionTimes(PeptideIdentificationList&, ...)` | `transform_peptide_identifications` |
| `transformRetentionTimes(IdentificationData&, ...)` | `transform_identification_data` |
| private `applyToFeature_`, `applyToConsensusFeature_` | `transform_feature`, `transform_consensus_feature` (public) |
| private `applyToBaseFeature_`, `storeOriginalRT_` | private plans in `alignment_transformer.rs` |
| the `store_original_rt` argument | the `store_original_rt` field |

`tests/alignment_transformer.rs` reuses the four container examples from the source transformer test and checks preservation of original RTs, chromatogram arrays, decreasing transforms, subordinate hulls, consensus handles, assigned/unassigned peptide RTs, the optional spectrum-ID extension, limits, and late-error atomicity; the identification-data overload, which the class test does not exercise, is checked with the peptide section's anchors and sequence (derived).
