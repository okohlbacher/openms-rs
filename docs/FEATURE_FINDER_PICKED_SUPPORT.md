# FeatureFinderAlgorithmPicked

This document covers the port of `FEATUREFINDER/FeatureFinderAlgorithmPicked.h`
and its implementation at core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
The source class is the centroided-peptide feature finder that
FeatureFinderCentroided runs.

The port was staged over two packages of the early TOPP bundle. B6-FFAP-SEEDS
ported the front half — the parameters, the input checks of `run`, the
intensity, trace and isotope-pattern scores, the isotope-pattern precalculation
and seed selection (steps 0 to 3.2 of `run_`). B7-FFAP-FEATURES ported the back
half — the isotope fit of each seed, the mass-trace extension, the trace fit,
the cropping, the quality checks, the parallel seed loop and the overlap
resolution (steps 3.3 and 4). `run` now produces features end to end.

| Rust file | Content |
| --- | --- |
| [`algorithm.rs`](../src/analysis/feature_finder_picked/algorithm.rs) | defaults, `Settings`, `validate_input`, `Limits`, `Options`, `run`, `feature_stage` (the seed loop) |
| [`scoring.rs`](../src/analysis/feature_finder_picked/scoring.rs) | `position_score`, `nearest_from`, `IntensityThresholds`, `ScoreArrays`, `find_isotope`, `isotope_score` |
| [`seeds.rs`](../src/analysis/feature_finder_picked/seeds.rs) | `IsotopeWindows`, `SeedStage`, `ChargeSeeds`, `UserSeed`, `overall_score` |
| [`extension.rs`](../src/analysis/feature_finder_picked/extension.rs) | `OverallScores`, `find_best_isotope_fit`, `extend_mass_traces`, `extend_mass_trace` |
| [`fitting.rs`](../src/analysis/feature_finder_picked/fitting.rs) | `FittedModel`, `crop_feature`, `check_feature_quality`, `build_feature`, the abort reasons |
| [`resolution.rs`](../src/analysis/feature_finder_picked/resolution.rs) | `intersection`, `resolve_overlaps`, `annotate_apex` |

The helper types come from
[`helper_structs.rs`](../src/analysis/feature_finder_picked/helper_structs.rs)
([support](FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md)). The two trace
fitters come from
[`gauss_trace_fitter.rs`](../src/analysis/feature_finder_picked/gauss_trace_fitter.rs)
and
[`egh_trace_fitter.rs`](../src/analysis/feature_finder_picked/egh_trace_fitter.rs)
over the shared
[`trace_fitter.rs`](../src/analysis/feature_finder_picked/trace_fitter.rs)
([support](TRACE_FITTER_SUPPORT.md)). The binary32 isotope patterns come from
`src/chemistry/isotopes.rs` ([support](ISOTOPE_SOURCE_PRECISION_SUPPORT.md)).
None of those files is edited here.

Tests: [`tests/feature_finder_picked_seeds.rs`](../tests/feature_finder_picked_seeds.rs)
(steps 0 to 3.2) and
[`tests/feature_finder_picked.rs`](../tests/feature_finder_picked.rs)
(steps 3.3 and 4).
Manifest: [`tests/data/feature_finder_picked_provenance.json`](../tests/data/feature_finder_picked_provenance.json).
Nothing is feature-gated in the library. Both integration tests need `mzml` and
`paramxml`; the user-seed cases also need `featurexml`.

**Module edges.** The two packages add three cross-module edges, all of which
the recorded graph already held, so `tools/check_module_cycles.py` reports no
new edge:
- `analysis -> math`, for `pearson_correlation_coefficient`;
- `analysis -> param`, for `Param` and `DefaultParamHandler`;
- `analysis -> concept`, for `PROTON_MASS_U`, `UserParam::NUM_OF_DATAPOINTS`
  and `parallel::{Threads, map_collect}`.

`analysis -> kernel` and `analysis -> metadata` were already there; B7 adds no
further edge.

## API mapping

Every member of the header is listed.

### `FeatureFinderAlgorithmPicked`

| Source member | Rust | Notes |
| --- | --- | --- |
| base `DefaultParamHandler` | `Settings::from_parameters` over `default_parameters` | the handler is built per call; no stored `param_` |
| base `ProgressLogger` | not ported | progress is off by default in the source; nothing is printed |
| `MapType`, `SpectrumType`, `FloatDataArrays` | `MSExperiment`, `MSSpectrum`, `ScoreArrays` | see *Score arrays* below |
| `PeakType`, `Seed`, `MassTrace`, `MassTraces`, `TheoreticalIsotopePattern`, `IsotopePattern` (protected) | `Peak1D` and the `helper_structs` types | |
| `FeatureFinderAlgorithmPicked()` | `default_parameters`, `HANDLER_NAME` | |
| `setSeeds(const FeatureMap&)` | `seeds` argument of `run` and `SeedStage::run` | |
| `setData_(MSExperiment&&, FeatureMap&)` (private) | `run` consumes the experiment and returns `RunOutput` | |
| `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap&)` | `run`, `run_with_options` | returns `RunOutput` with the features, the log and the abort counts |
| `getDefaultParameters()` | `default_parameters()` | 29 entries, 7 section descriptions |
| `run_()` (protected) | `SeedStage::compute` (steps 0 to 3.2) then `feature_stage` (steps 3.3 and 4) | |
| `map_` | `SeedStage::experiment`, `into_experiment` | |
| `features_` | `RunOutput::features` | |
| `log_`, `debug_` | `Settings::write_debug` | debug output refused (native differences) |
| `aborts_`, `abort_()` | `RunOutput::aborts`, aggregated serially | the counts are exact; the source races them (candidate 5) |
| `abort_reasons_` | not ported | debug only, and the debug output is refused |
| `seeds_` | `SeedStage::user_seeds` (`UserSeed`: m/z and RT) | the only fields the source reads |
| `pattern_tolerance_` ... `max_feature_intersection_`, `reported_mz_` | `Settings` fields of the same names; `reported_mz` as `ReportedMz` | |
| `intensity_rt_step_`, `intensity_mz_step_`, `intensity_thresholds_` | `IntensityThresholds::rt_step`, `mz_step`, `quantiles` | |
| `isotope_distributions_` | `IsotopeWindows::patterns` | |
| `updateMembers_()` | `Settings::from_parameters` | also reads the `run_` locals: charges, `fit:max_iterations`, abundances, seed thresholds, user-seed tolerances, `feature:min_score`, `feature:rt_shape` |
| `intersection_()` | `resolution::intersection` | |
| `getIsotopeDistribution_(double)` | `IsotopeWindows::get` | |
| `findBestIsotopeFit_()` | `extension::find_best_isotope_fit` | returns the score and the pattern instead of an out-parameter |
| `extendMassTraces_()`, `extendMassTrace_()` | `extension::extend_mass_traces`, `extension::extend_mass_trace` | the overall-score array is `extension::OverallScores` |
| `nearest_()` | `nearest_from` | also returns the steps walked |
| `findIsotope_()` | `find_isotope` | returns work units; checks indices |
| `positionScore_()` | `position_score` | |
| `isotopeScore_()` | `isotope_score` | |
| `intensityScore_(Size spectrum, Size peak)` | `IntensityThresholds::score(rt, mz, intensity)` | |
| `intensityScore_(Size rt_bin, Size mz_bin, double)` | `IntensityThresholds::bin_score` | `None` outside the grid |
| `chooseTraceFitter_(double&)` | `fitting::FittedModel::new` from `Settings::rt_shape` | the enum replaces the pointer plus the `tau != 0` flag and the `dynamic_pointer_cast` |
| `cropFeature_()` | `fitting::crop_feature` | returns the cropped traces instead of an out-parameter |
| `checkFeatureQuality_()` | `fitting::check_feature_quality` | returns `QualityOutcome`: `Accepted(FeatureQuality)` or `Rejected(reason)` |
| step 3.3 (`.cpp:576-856`), the `omp parallel for` and the containment pass | `feature_stage`, with `fitting::build_feature` for step 3.3.5 | `concept::parallel::map_collect` over the seed indices |
| step 4 (`.cpp:859-1016`) | `resolution::resolve_overlaps`, `FeatureMap::sort_by_mz`, `retain`, `sort_by_intensity(true)`, `resolution::annotate_apex` | |
| `writeFeatureDebugInfo_()` | not ported | debug output refused |
| `operator=`, copy constructor (private, not implemented) | not applicable | the stage is an owned value |
| step 1, second half (`.cpp:280-287`) | `fill_intensity_scores` (crate-private) | |
| step 2 (`.cpp:291-348`) | `fill_trace_scores` (crate-private) | |
| step 3.1 (`.cpp:447-488`) | `fill_pattern_scores` (private) | |
| step 3.2 (`.cpp:489-574`) | `select_seeds` (private), `overall_score`, `SeedStage::charges` | |

Native additions: `Limits`, `Options`, `AbundanceOverride`, `RunOutput`,
`validate_input`, `UNSORTED_WARNING`, `BASE_MAX_ISOTOPES`,
`OVERRIDE_EXTRA_ISOTOPES`, `Settings::charge_count`, `Settings::max_isotopes`,
`QUANTILE_COUNT`, `ScoreArrays`, `ChargeSeeds`, `UserSeed`,
`IsotopeWindows::mass_window_width`, the `SeedStage` accessors, and from the
back half `OverallScores`, `FittedModel`, `FeatureQuality`, `QualityOutcome`,
`FeatureInput`, the seven `ABORT_*` reason constants, `SPECTRUM_INDEX`,
`SPECTRUM_NATIVE_ID` and `invalid_apex_warning`.

### `FeatureFinderDefs` (same header)

`IndexPair`, `ChargedIndexSet`, `IndexSet`, `enum Flag { UNUSED, USED }` and the
exception class `NoSuccessor` are not ported. The algorithm uses none of them.
The struct duplicates `FEATUREFINDER/FeatureFinderDefs.h`
([helper structs support](FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md)).

## Preserved source conventions

- **Parameters.** The defaults carry the source names, values, descriptions,
  numeric and string restrictions and `advanced` tags. Their insertion order is
  the source's, and they equal the executed `-write_ini` section. Unknown
  entries are warnings, and type or restriction violations are errors, as in
  `DefaultParamHandler::setParameters`. `min_spectra` is
  `floor(mass_trace:min_spectra * 0.5)`. The three percentages are divided by
  100. An abundance is changed when its `ParamValue` differs exactly from the
  default. Each changed abundance adds 1000 isotopes to the 20 of the patterns.
- **Check order of `run`.** The checks run in the source's order:
  1. no spectra returns an empty map before the parameters are read;
  2. no peak in spectra or chromatograms is an error;
  3. MS levels other than `{1}` are an error;
  4. unsorted spectra are sorted, with the source warning;
  5. a negative first m/z is an error.

  The messages are the source's. The parameters are applied after these checks.
- **Intensity bins (step 1).** The grid spans the MS1 RT and m/z ranges.
  - Bin borders are `start + i * step` and `start + (i + 1) * step`, evaluated
    in that form, and both are inclusive, as `areaBeginConst` is.
  - Intensities are promoted to `f64` and sorted. Quantile `i` is element
    `floor(0.05 * i * (n - 1))`, and an empty cell keeps 21 zeros.
  - A peak's half-bin position is `floor((x - start) / step * 2)`, capped at
    `2 * bins - 1`, and selects the neighbouring bins with the source's
    edge, odd and even rules.
  - The four cell scores are weighted by `sqrt((1 - d_rt)^2 + (1 - d_mz)^2)`
    over their sum, in the source's order.
  - The cell score uses `lower_bound` and the `0.05` interpolation. The clamp
    to `[0, 1]` lets NaN through.
  - Scores are stored as `f32`.
- **Trace scores (step 2).** The scans `min_spectra..len - min_spectra` are
  scored; the others keep zero.
  - The following scans are visited before the preceding ones, empty scans are
    skipped, and the nearest peak is found with the lower-m/z midpoint rule.
  - Position scores are summed in `f64` and divided by `2 * min_spectra`.
  - A peak stops being a local maximum when a neighbour with a positive
    position score is strictly more intense in `f32`.
- **Isotope windows (step 2.5).**
  - There are `ceil(max_mz * charge_high / width) + 1` windows. Window `i`
    estimates the peptide mass `0.5 * width + i * width`, with at most 20 (or
    1020, 2020) isotopes, in binary32 source precision.
  - The source `trimLeft` keeps a pattern whole when no weight reaches the
    cutoff. It is followed by `trimRight`, both at
    `intensity_percentage_optional`; `trimmed_left` counts the removed leading
    isotopes.
  - Optional isotopes are counted with the source's `is_begin`/`is_end` state
    machine against `intensity_percentage`.
  - The maximum is taken with a strict `>` from 0, and every intensity is
    divided by it.
- **Pattern scores (step 3.1).**
  - The pattern is anchored at the first maximum (`std::max_element`), with
    positions `mz + (i - max_isotope) / c`. The walk starts at the peak nearest
    to `mz - (len + 1) / c`.
  - `nearest_` walks upwards while strictly closer.
  - `findIsotope_` matches on a non-zero position score, so NaN matches, in the
    centre, previous and next scan in that order. The found peak comes from the
    centre, or else from the first neighbour that matched, and the intensities
    and position scores are averaged.
  - `isotopeScore_`:
    - a missing required peak scores 0;
    - the search starts behind the last missing optional peak at each end;
    - two isotopes are allowed only for the starting candidate;
    - NaN correlations count as 0, two-isotope fits are capped at
      `min_isotope_fit`, and a candidate must improve by `1 +
      optional_fit_improvement` over a best that starts at 0.01;
    - the best trailing count is re-read at the start of every inner loop, so a
      new best fit narrows the later candidates;
    - left-out peaks become `Removed`, and the m/z factor is the mean m/z score
      of the kept isotopes.
  - A peak's pattern score is raised only by a higher score (compared in `f64`,
    stored as `f32`) and only for matched, non-removed peaks.
- **Seeds (step 3.2).**
  - The overall score is the `f32` product `trace * intensity * pattern`,
    formed left to right, raised to `1.0f / 3.0f`. It is stored for every peak
    of the scored scans and compared with the thresholds in `f64`.
  - Seeds are local maxima at or above `seed:min_score`. With user seeds, the
    threshold is `user-seed:min_score`, and a peak also needs a user seed
    strictly within both user-seed tolerances. The search is the source's
    `lower_bound` plus forward walk over the seeds sorted by m/z.
  - Seeds are sorted by descending `f32` intensity with
    `Seed::is_less_intense_than`.
  - Charges are processed from `charge_low` upwards, and the log line is
    `Found <n> seeds for charge <c>.`.
- **Isotope fit (step 3.3.1, `findBestIsotopeFit_`).**
  - The theoretical pattern is the window of `seed_mz * charge`, and the search
    window is `(isotopes + 1) / charge` wide on either side of the seed peak,
    found by the source's two linear walks with their off-by-one corrections
    (`--end`, `++begin`).
  - Every peak of that window starts a placement whose isotope `i` sits at
    `start_mz + i / charge`, matched with `findIsotope_` from a peak index that
    the previous isotope left behind.
  - A placement counts only while it still matches the seed's own peak, checked
    before *and* after `isotopeScore_`, which may drop optional isotopes.
  - The score is `isotopeScore_` **without** the m/z-distance factor, and the
    best is kept with a strict `>`, so the earliest of equal scores wins.
  - `best_pattern.theoretical_pattern` is assigned after the loop whether or not
    a placement qualified.
- **Mass-trace extension (step 3.3.1, `extendMassTraces_`).**
  - The maximum trace starts at the strongest matched isotope, compared as
    `float > double` against a running maximum that starts at 0, so a
    zero-intensity peak never starts it.
  - It is extended without bounds; its first and last retention time then bound
    every other trace. The candidate is dropped when it has fewer than three
    peaks or fewer than `2 * min_spectra - max_missing` peaks, the latter in
    wrapping `UInt` arithmetic.
  - Every other matched isotope first looks for a stronger start peak in the
    scans `spectrum - min_spectra .. spectrum + min_spectra`, where the lower
    bound is an unsigned difference that wraps for the first `min_spectra`
    scans and empties the range. Inside that search the nearest peak is looked
    up at the *current* start peak's m/z, which moves, while the tolerance is
    compared against the *original* m/z (candidate 2).
  - An isotope removed or missing during the isotope fit is skipped.
  - A trace with fewer than three peaks is handled by comparing the **pattern**
    index with `MassTraces::max_trace` (candidate 1).
- **One trace (`extendMassTrace_`).**
  - The downward call walks from `spectrum - 1`; the upward call reverses the
    peaks collected so far and walks from `spectrum + 1`, so the trace ends in
    chronological order.
  - A peak is missing when none was found, its overall score is below `0.01`, or
    its m/z score is exactly zero; every comparison is false for NaN, so a NaN
    counts as present. More than `max_missing` consecutive missing scans stop
    the walk.
  - The slope bound is doubled when retention-time bounds are given
    (`max_rt != min_rt`, which is how the source detects them).
  - The intensity deltas start as `min_spectra - 1` zeros; each accepted peak
    appends `(intensity - last) / last` in `f64`, and the mean of the last
    `min_spectra` deltas, summed left to right from `0.0`, is compared with the
    bound. Exceeding it removes the last `min(added, min_spectra - 1)` peaks of
    *this* call and stops the walk.
- **Fit, cropping and checks (steps 3.3.2 to 3.3.4).**
  - The baseline is `0.75 *` the lowest peak intensity over all traces, and the
    maximum trace's maximum peak is refreshed, both before the fit.
  - `feature:rt_shape` picks the Gaussian or the EGH model, which receives only
    `max_iteration` from `fit:max_iterations`, so the tool path always fits
    unweighted.
  - Cropping keeps the peaks inside `[lower_rt_bound, upper_rt_bound]`
    inclusive, scores each trace by
    `sqrt(max(0, correlation) * max(0, 1 - mean relative deviation))` and drops
    bad traces by the source's position rules relative to `max_trace`. The
    baseline is copied to the result last.
  - The five quality rules run in the source's order: `checkMaximalRTSpan`,
    `MassTraces::isValid`, the centre inside the retention-time bounds,
    `checkMinimalRTSpan`, and `sqrt(correlation * fit_score)` against
    `feature:min_score`. Each `std::max(0.0, x)` is the source's
    `(0.0 < x) ? x : 0.0`, so a NaN becomes zero.
- **Feature creation (step 3.3.5).**
  - Fields in the source's order: the label (`MetaInfoRegistry` index 3, the key
    `label`), the charge, the overall quality (narrowed to `f32`), `score_fit`,
    `score_correlation`, the retention time (the model's centre), the width (the
    model's FWHM, narrowed to `f32`, which also writes the `FWHM` meta value),
    `num_of_datapoints`, the three `EGH_*` values for an asymmetric fit, the
    m/z, the intensity, and one convex hull per trace.
  - `feature:reported_mz` selects the maximum trace's average m/z, the
    intensity-weighted average over the cropped traces, or the monoisotopic
    value `average - PROTON_MASS_U / charge * (max trace index + trimmed_left)`
    (candidate 3).
  - The intensity is the model's area divided by the maximum of the isotope
    window of the feature's **m/z**, not of its mass (candidate 4).
  - The seeds swallowed by the feature are the later seeds inside both the
    overall bounding box (inclusive) and one of the mass-trace hulls.
- **The seed loop and the containment pass (step 3.3).**
  - `plot_nr` counts the seeds that reached the fit; the label of every accepted
    feature is then overwritten with the running feature number, which continues
    across charges.
  - The containment pass runs serially in seed order: a candidate whose seed a
    previously accepted feature swallowed is dropped, and the seeds of an
    accepted feature are added to the swallowed set.
  - The log line is `Found <n> feature candidates for charge <c>.`, printed
    directly after that charge's seed line.
- **Overlap resolution (step 4).**
  - The map is sorted by m/z; each feature's overall bounding box and the
    largest m/z extent are precomputed once.
  - The inner loop stops at `f2.mz - f1.mz > 2 * max_mz_span`, skips pairs with
    a zero intensity or disjoint boxes, and acts when
    `intersection_ >= feature:max_intersection`.
  - Equal charges compare the `f32` product `intensity * quality`; otherwise a
    charge that is a multiple of the other keeps the higher charge; otherwise
    the higher quality wins. A tie keeps the first feature. The loser is cloned
    into the winner's subordinates with the intensity it had, and its own
    intensity becomes zero.
  - Zero-intensity features are then removed, the map is sorted by descending
    intensity, and each feature receives `spectrum_index` (`RTBegin(rt)`) and,
    when that index addresses a scan, `spectrum_native_id`.

## Native differences

| Source behaviour | This port | Reason and evidence |
| --- | --- | --- |
| scores are float data arrays appended to each spectrum, replacing its existing float arrays | `ScoreArrays`: one flat `f32` array per score, outside the spectra | Only the algorithm reads them, and only debug mode writes them. The input spectra keep their arrays, and no per-spectrum allocation is needed. |
| `spectrumRanges().byMSLevel(1)` needs `updateRanges()` from the caller | ranges are computed on demand from the validated spectra | No stale-range state exists, so the source's FAIMS "No ranges for this MS level" crash cannot occur. The source message "needs updated ranges" belongs to the peak-count check and is kept verbatim. |
| NaN or infinite RT, m/z or intensity values are sorted and binned with undefined results | `Error::InvalidValue` | The native readers never produce them. |
| `mass_trace:min_spectra = 1` gives `min_spectra_ = 0`. Every trace score becomes 0/0 = NaN and every peak a local maximum; no overall score reaches a threshold, no seed is found and the run returns an empty map | the same: NaN trace scores, no seed, an empty map | The execution (B6 driver, `ffc1_min_spectra_1`) shows the source is *defined* here, so the port follows it (lead decision of 2026-09-15, `CPP-271`). B6 refused the configuration; that refusal is gone. Nothing later in the algorithm is reached, so the source's `size_t(-1)` delta buffer in `extendMassTrace_` stays unreachable; the port returns `Error::InvalidValue` if it ever is. |
| a single retention time or a single m/z makes the bin width 0; `(UInt)floor(NaN)` is undefined | `Error::InvalidValue` | Undefined behaviour. |
| `charge_low > charge_high + 1` wraps the `UInt` charge count and indexes past the score arrays | `Error::InvalidValue` from `Settings::charge_count` | Undefined behaviour. `charge_low == charge_high + 1` gives zero charges, as in the source. |
| `write_debug = true` writes `debug/` into the working directory and then throws on the undeclared `debug:pseudo_rt_shift` | `Error::Unsupported` | A defect, and library code does not write files. |
| a changed `abundance_12C` or `abundance_14N` builds the override from a default `IsotopeDistribution` that already holds `(0, 1)`; the patterns grow (FFC_1 window 0: 27 normalised bins) and FFC_1 with 12C = 90 % finds 0 seeds, 0 candidates and 0 features (C2 `ffap_ffc1_abundance_12C_90`) | the intended two-isotope distribution, which **does** find seeds and features | `CPP-247`, lead decision of 2026-09-15: follow the intent, not the defect, because the generator rejects the stray-peak construction and refusing a parameter the source accepts is worse. **This is the one place where the port's features differ from the executed C++ by design.** `AbundanceOverride::Refuse` is the opt-in for a caller that must not diverge. |
| the overall score is `std::pow(float, float)`, the platform `powf` | `libm::pow` in `f64`, rounded once to `f32` | Apple `powf` misrounds 99 of the 30,840 retained overall scores by one binary32 step, 12 of 3,084 on FFC_1. The port's value is the correctly rounded power for all of them (60-digit decimal check) and is the same on every machine. `libm::powf` differed on 226 of 3,084. No seed list changes. `tests/data/feature_finder_picked/overall_rounding.tsv` lists every difference. |
| `std::sort` of seeds with equal `f32` intensity is unspecified | stable: scan, then peak order | Deterministic. No retained configuration has a tie (the drivers check adjacent ties). |
| `std::sort` of a bin's intensities: `-0.0` and `+0.0` compare equal | `total_cmp`: `-0.0` first | Only a quantile's zero sign can differ, and only for inputs with both signed zeros; the tool path filters non-positive intensities. |
| progress, `Found N seeds for charge c.` and `Found N feature candidates for charge c.` go to `std::cout`; the overlap count, the abort reasons, the feature count and the apex warning to `OPENMS_LOG_INFO`/`WARN` | `SeedStage::log` and `RunOutput::log`, in the source's order | Library code never prints. The candidate line is inserted directly after its charge's seed line, so the two `std::cout` lines are adjacent as in the source, even though this port computes every charge's seeds first. The bare newlines the source logs around the abort block (`FeatureFinderAlgorithmPicked.cpp:1019` and `1026`) are emitted as empty log entries, so a caller that prints the log line by line — `FeatureFinderCentroided` does — reproduces the executed C++ stdout block exactly. |
| steps 3.1 to 3.3 run per charge | steps 3.1 and 3.2 run for every charge first | Step 3.3 reads only its own charge's arrays, so the arrays and seeds are identical. |
| user seeds are a copied `FeatureMap` sorted with `std::sort` | positions only, sorted stably; non-finite positions are `Error::InvalidValue` | NaN breaks the source sort. Equal m/z values give the same search result in any order. |
| unbounded work | `Limits`: spectra, peaks, charges, bins per dimension, windows, pattern values, score bytes, work units | Checked before the allocation or computation each bounds. The FFC_1 workload is several orders of magnitude below every default. |
| `Math::pearsonCorrelationCoefficient` returns an infinity when its denominator underflows to 0 from non-zero deviations | NaN, which counts as 0 | Inherited from `src/math/statistic_functions.rs`. Unreachable at isotope intensity scales. |
| `Exception::UnableToFit` from `fitter->fit` is thrown inside the `omp parallel for`, where it is not caught, so the whole run ends | the failure becomes that seed's abort reason and the run continues | The source's own handling of every other failure at this point is an abort, and an uncaught exception leaving an OpenMP region is undefined. No executed configuration reaches it. |
| `extendMassTraces_` dereferences `pattern.spectrum[0]` when the pattern matched no peak | `Error::InvalidValue` | Undefined behaviour. Unreachable from `run_`, where the pattern always contains the seed; reachable through the public function, which the test exercises. |
| `traces[traces.max_trace]` is indexed before the fit without a range check | `Error::InvalidValue` | Undefined behaviour when `max_trace` is stale. The one branch that could make it stale (`traces.clear()` for a trace before `max_trace`) is unreachable, because `max_trace` is still 0 at every index that could satisfy `p < max_trace`. |
| `setWidth` stores any FWHM, including a NaN produced by a non-finite fit that the quality checks let through (every comparison against a NaN is false) | `Error::InvalidValue` from `BaseFeature::set_width` | The kernel setter validates. No executed configuration produces a non-finite fit. |
| `f2.getCharge() % f1.getCharge()` divides by zero for a zero charge | the branch is skipped and the quality rule decides | `isotopic_pattern:charge_low` has a minimum of 1, so a zero charge cannot reach step 4; a trap would be worse than the fall-through. |
| a hull with no point has the default `DBoundingBox` `[DBL_MAX, -DBL_MAX]`, whose `width()` is negative infinity and poisons `intersection_` | such a hull is skipped | Every hull built here holds at least three points. |
| `plot_nr` is assigned in an OpenMP critical section, so its value depends on the schedule | assigned in seed order | It is overwritten by the feature number for every feature that survives, and only the refused debug output reads it otherwise. |
| `aborts_[reason]++` runs inside the parallel region without synchronisation | aggregated serially in seed order | A data race (candidate 5). The port's counts are exact and independent of the thread count; the C2 driver records the library's map single-threaded only for the same reason. |
| the seed loop is an OpenMP `parallel for` with four named critical sections | `concept::parallel::map_collect` over the seed indices with `Options::threads` | `map_collect` returns results in input order, so no critical section is needed and the output is bit-identical at 1, 2 and 8 threads (`the_seed_loop_is_bit_identical_across_thread_counts`). The source's results are schedule-independent for the same reason, its `tmp_feature_map` being a `std::map` keyed by seed index; the C++ oracle gave identical output at 1, 2 and 8 threads. |
| the containment pass scans a growing `std::vector<Size>` of swallowed seeds | a `BTreeSet` | A pure membership test; duplicates in the source's vector change nothing. |
| `FeatureMap::sortByMZ` and `sortByIntensity` use `std::sort`, which is unstable | stable sorts | Only exactly tied m/z or intensity values could be ordered differently. None of the six executed configurations has one. |
| `setMetaValue("spectrum_index", Size)` stores an unsigned value | `i64`, refused above `i64::MAX` | The source's `DataValue` narrows to a signed integer anyway; the featureXML writer writes the same digits. |
| `CoarseIsotopePatternGenerator` iterates elements in heap-address order, which varies between runs | ascending atomic number | Inherited from B2. It is the majority order (198 of 200 runs), which every retained execution used; all windows match. |
| `MSExperiment::sortSpectra` | `MSExperiment::sort_spectra` validates each spectrum first | Existing kernel API; a spectrum with inconsistent attached records is an error. |

### Score arrays

The source names are kept by `ScoreArrays::array_names`: `trace_score`,
`intensity_score`, `local_max`, `pattern_score_<c>` for each charge, then
`overall_score_<c>` for each charge. Trace scores, local-maximum flags and overall
scores stay zero outside `min_spectra..len - min_spectra`, as in the source.

## Checked boundaries and evidence

### Tier 1: executed C++ (product SDK, Debug, core `4fdec46`)

The traced sources are hash-identical between `4fdec46` and the pin. Every float
of the **seed stage** is compared bit for bit on Linux x86-64; the only exception
there is the documented overall-score rounding, whose correctly rounded values
the test asserts instead. The **feature stage** compares counts, identities and
orders exactly and fitted quantities within `1e-9` relative, because the
Levenberg-Marquardt transcription is not yet bit-faithful to the executed Eigen;
the measured agreement is in the next section.

| Evidence | Configurations | Test |
| --- | --- | --- |
| C1 `FFC_write_ini`: FeatureFinderCentroided `-write_ini` (two runs, identical) | the whole algorithm section | `default_parameters_equal_the_executed_write_ini_section` |
| C2 `ffap_stages` library state (`omp1` and `omp4`, identical): effective members, bins, 21 quantiles, 15 windows, all score arrays of 3,084 peaks, printed seed counts | FFC_1 INI (tool loading), class-test INI, #9247 tight pattern and tight trace, FFC_1 with the retained output as user seeds | `ffc_1_stage_...`, `class_test_stage_...`, `tolerance_swap_stages_...`, `user_seed_stage_matches_the_executed_library` |
| B6 `seed_stage` library state (two runs, byte-identical), same content | default parameters (10 bins, charges 1 to 4, width 25, 107 windows); FFC_1 INI with 7 bins and charges 1 to 3 (21 windows) | `default_parameter_stage_...`, `seven_bin_three_charge_stage_matches_the_executed_library` |
| B6 `seed_stage`, `ffc1_min_spectra_1` | `mass_trace:min_spectra = 1`: `min_spectra_` 0, 0 seeds, 0 features, exit 0 | `min_spectra_one_follows_the_source_and_finds_no_seed` |
| C2 `ffap_stages` **final `FeatureMap`**, `aborts_` and the two `std::cout` lines, for six configurations | FFC_1 symmetric, FFC_1 asymmetric (EGH), FFC_1 with the retained output as user seeds, the class-test INI, the two `#9247` tolerance swaps | `every_configuration_matches_the_executed_library` |

The FFC_1 score table also pins every loaded peak's m/z and intensity bits
against the C++ loader.

The oracle has two single-bin configurations (FFC_1 and the class test). The
B6 driver adds 7 and 10 bins, which exercise the four-cell interpolation, and up
to four charges.

### The feature stage: what was compared and what it showed

`tests/feature_finder_picked.rs` reproduces each of the six executed
configurations twice: once through `run`, against the library's final
`FeatureMap` (tier 1), and once step by step through the public functions of
`extension.rs` and `fitting.rs`, against the C2 driver's per-seed replay
(adapted). 4,649 numeric values are compared in total.

Measured on both platforms, Linux x86-64 (IBMI dax, glibc) and macOS arm64. The
two columns agree everywhere except inside the two fits named below, so the
counts are given once.

| Family | Values | Bit-identical | Largest relative departure |
| --- | ---: | ---: | --- |
| isotope fit score (`findBestIsotopeFit_`) | 131 | 131 | 0 |
| isotope pattern intensities and m/z scores | 1,572 | 1,572 | 0 |
| mass traces: peak identities, `theoretical_int`, baseline | 1,629 | 1,629 | 0 |
| fitted model (centre, height, FWHM, area, bounds, sigma, tau) | 739 | 176 | `2.3e-3`, two seeds only (see below); `6.5e-10` everywhere else |
| quality scores (`fit_score`, `correlation`, `final_score`) | 291 | 51 | `3.2e-10` |
| final feature fields | 287 | 206 | `3.2e-10` |

Peak identities are `(spectrum, peak)` index pairs and are compared as text, so
"bit-identical" there means the same peaks in the same order. Counts, charges,
labels, `num_of_datapoints`, hull counts, hull point counts, hull point
coordinates, subordinate counts, abort reasons and abort counts are compared
exactly and agree everywhere. Of the final feature fields, every m/z, intensity,
quality, width and `FWHM` of all 33 features across the six configurations is
bit-identical; the retention time departs by at most `7.7e-13`,
`score_correlation` by at most `1.1e-11`, `score_fit` by at most `3.2e-10` and
the three `EGH_*` values by at most `5.3e-16`.

**The one departure beyond `1e-9`.** Seeds 11 and 12 of
`classtest_9247_tight_pattern` fit the same traces, and their fitted parameters
depart from the executed Eigen by up to `2.3e-3` on Linux (the area; sigma and
FWHM `1.3e-3`, height `9.3e-4`, upper bound `6.7e-5`, centre `2.6e-5`, lower
bound `1.7e-5`) and by up to `5.8e-4` on macOS arm64 (area; sigma and FWHM
`3.4e-4`, height `2.4e-4`, centre `7.3e-6`). These two fits are also the only
results that differ between the two platforms; every other compared value agrees
on both to within `6.5e-10`. The cause is the port's
Levenberg-Marquardt transcription, not this package: the *inputs* of those two
fits — the peak identities, the theoretical intensities and the baseline — are
bit-identical to the executed ones, and only the solver's output differs. It is
the gap `docs/TRACE_FITTER_SUPPORT.md` records under "Known gap: solver fidelity
beyond the fixtures", whose agreement was measured on the GaussTraceFitter and
EGHTraceFitter class tests and on the FeatureFinderCentroided_1 seeds, neither
of which covers this configuration. Lane B3b is root-causing it. Both seeds are
rejected by `checkFeatureQuality_` in the executed C++ *and* here, with the same
reason, so no feature changes. The test pins the larger measured bound in
`KNOWN_FIT_GAP` and asserts that such a seed never becomes a feature; the
tolerance for every other seed stays `1e-9`.

### Adapted

The ordered seed lists come from replicas of `.cpp:493-548`, in C2 `ffap_stages`
and B6 `seed_stage`. The replicas apply the source selection code verbatim to
the library's score arrays. C2 checks that recomputing the overall scores
reproduces the library arrays with 0 mismatches. The seed counts of both
replicas equal the library's printed `Found N seeds` lines, which the tests
assert against the log. The lists are exact (spectrum, peak, `f32` intensity,
overall score) for all 7 configurations: 25, 25, 15, 18 and 24 seeds for charge
2, 42/44/0/0 for charges 1 to 4 and 16/23/0 for charges 1 to 3.

### Tier 3 and 4

- **FeatureFinderAlgorithmPicked_test.** All four sections are covered. The
  constructor and destructor sections by `constructor_and_destructor`; the class
  test has no `getDefaultParameters` section of its own, and the executed
  `-write_ini` comparison covers it. The `run()` section's literals — 8
  features, `num_of_datapoints` 88/71/47, the eight qualities at absolute 0.001
  and the eight intensities at absolute 20.0 — are asserted by
  `class_test_run_finds_the_expected_eight_features`, and the `[EXTRA] #9247`
  section's 1-and-0 features with the surviving feature's RT, m/z, quality and
  intensity by `class_test_tolerance_swap_is_directional`.
- **Source review.** Covered cases:
  - the input check order, including invalid parameters on empty input and on
    MS2 input;
  - unsorted input: sorted with the warning, and scored exactly like the sorted
    input;
  - the refusals: `write_debug`, charge wrap, zero-width ranges, non-finite
    user seeds, abundances;
  - the conversions in `settings_follow_update_members`;
  - restriction and type violations.
- **Hand-derived (feature stage).**
  - `intersection_`: the two containment cases, both partial-overlap cases,
    disjoint and touching boxes, the division by the smaller total, and that a
    hull pair whose boxes do not intersect contributes nothing
    (`intersection_follows_the_source_cases`);
  - an isotope pattern that matched no peak is refused instead of dereferenced
    (`an_empty_pattern_is_refused_instead_of_dereferenced`);
  - the cropping position rules on a two-trace candidate whose second trace lies
    beyond the model's bounds (`cropping_follows_the_source_position_rules`);
  - the seed-loop ceilings fail below the FFC_1 workload and pass at the
    defaults (`the_seed_loop_ceilings_are_checked_first`);
  - the labels are `0..n` after the containment pass and every feature carries
    both apex meta values (`labels_are_the_feature_numbers_in_order`).
- **Hand-derived (seed stage).**
  - `position_score` at both branches and with a zero tolerance;
  - `nearest_from` walks, ties and the first local minimum;
  - `find_isotope` centre and neighbour matching and index errors;
  - `isotope_score`: required peaks, the two-isotope cap, the 0.01 start and
    the narrowed inner loop, on a pattern where the skipped candidate would
    have won;
  - `bin_score` interpolation and clamping;
  - `score` at the grid centre.
- **Area-iterator equivalence.** `intensity_bins_equal_the_area_iterator`
  compares every cell's quantiles for 1, 7 and 10 bins with the kernel's
  `MSExperiment::area_begin`.
- **Limits.** Each seed-stage ceiling fails one below the FFC_1 requirement and
  passes at it (`limits_are_checked_before_the_work`); the two seed-loop
  ceilings likewise (`the_seed_loop_ceilings_are_checked_first`).
- **Determinism.** `the_seed_loop_is_bit_identical_across_thread_counts` runs
  FFC_1 at 1, 2 and 8 threads and compares every feature's bits, metadata,
  hulls and subordinates, plus the log and the abort counts.

### Performance

FeatureFinderCentroided is the bundle's second benchmark target. The input is
FeatureFinderCentroided_1: 112 MS1 scans, 3,084 peaks, one charge, 25 seeds,
8 features.

| Measurement | Scope | Best of n |
| --- | --- | --- |
| Rust `run`, 1 thread (Linux, IBMI dax, AMD EPYC 9654, release) | the algorithm only | **9.6 ms** (9) |
| Rust `run`, 2 threads / 8 threads (same host) | the algorithm only | 7.4 ms / 5.5 ms |
| Rust mzML load (same host) | `FileHandler::load_experiment_with_options`, MS1 only | 2.4 ms (9) |
| C++ `FeatureFinderCentroided`, Release `bc9cc12-c19e494-174b576`, `-threads 1` and `OMP_NUM_THREADS=1` (same host) | TOPPBase startup, INI, mzML load, the algorithm and the featureXML write | 90 ms wall, 48 MB peak (7) |
| C++ `FileInfo` on the same mzML (same host) | TOPPBase startup, mzML load and the report | 40 to 70 ms, 42 MB peak (3) |
| Rust `run`, 1 thread (macOS arm64, release) | the algorithm only | 6.1 ms (7) |

The two tools are not directly comparable — the C++ number is end to end and
this port's wrapper is C5's and B10's — but subtracting the `FileInfo` load cost
leaves roughly 30 to 50 ms for the C++ algorithm against this port's 9.6 ms at
one thread. That subtraction is an estimate, not a measurement; the instrumented
tool-level comparison belongs to B10.

The source's algorithmic complexity is kept everywhere. The seed loop allocates
one `MassTraces` and one `IsotopePattern` per seed and reuses the pattern buffer
across the placements inside `findBestIsotopeFit_`, so nothing is allocated per
peak; the per-peak score arrays are three flat `f32` vectors plus two per charge,
allocated once. `map_collect` builds one rayon pool per charge, which is the only
threading cost.

### Verification commands

Run from the worktree root on the IBMI build node through
`~/.local/bin/openms-kim-gate.sh`; `target_verification` in the manifest records
the results.

```text
cargo fmt --all -- --check
cargo test --locked --all-features --test feature_finder_picked \
  --test feature_finder_picked_seeds --test feature_finder_picked_helper_structs \
  --test trace_fitter --test gauss_trace_fitter --test egh_trace_fitter \
  --test isotopes_source_precision --test mass_trace --test mass_trace_detection
cargo +1.85.0 test (same targets)
cargo test --locked --no-default-features --features mzml,paramxml,featurexml --test feature_finder_picked
cargo clippy --locked --all-features --all-targets -- -D warnings
cargo doc --locked --all-features --no-deps
python3 tools/check_core_sdk.py
python3 tools/check_module_cycles.py
python3 tools/check_doc_coverage.py
```

## C++ issue candidates

### From the feature stage (B7)

New candidates for the integrator's `OpenMS_CPP_ISSUES.md`; the numbering below
is local to this document, because the log's identifiers are the integrator's to
assign. All five are reproduced by this port, because they change ordinary
output.

1. **`extendMassTraces_` compares a pattern index with a trace index.** At
   `.cpp:1467-1475` an isotope whose trace stays below three peaks is handled by
   `p < traces.max_trace` / `p > traces.max_trace`, where `p` indexes the
   *isotope pattern* and `max_trace` indexes the *traces collected so far* and is
   still 0 until the maximum trace is reached. An invalid trace at `p == 0` is
   therefore appended instead of skipped, any later invalid trace before the
   maximum stops the extension, and the `traces.clear()` branch the comment
   describes ("Missing traces in the middle of a pattern are not acceptable") is
   unreachable: `p < max_trace` requires `p < max_trace_index`, where
   `max_trace` is still 0.
2. **The better-seed search reads a moving m/z.** At `.cpp:1425-1445` the
   nearest-peak lookup uses `map_[starting_peak.spectrum][starting_peak.peak]
   .getMZ()`, which is re-read after `starting_peak` moves, while the tolerance
   test compares against `mz`, captured before the loop. The search therefore
   drifts with the accepted candidates while the acceptance window does not.
3. **The monoisotopic m/z uses the proton mass and a trace index.** At
   `.cpp:782-784` the correction is `PROTON_MASS_U / charge * (trace index +
   trimmed_left)`. The spacing between isotope peaks is the neutron-ish
   `C13C12_MASSDIFF_U` (1.00335 u), not `PROTON_MASS_U` (1.00728 u), and the
   trace index is a position in the pattern rather than an isotope number. The
   reported monoisotopic m/z is biased by about 4 mDa per isotope step.
4. **The final intensity picks the isotope window by m/z, not by mass.** At
   `.cpp:790` the divisor is `getIsotopeDistribution_(f.getMZ()).max`, while the
   windows were precalculated over *mass* (`.cpp:356`, `max_mz * charge_high`).
   For a doubly charged feature the window is the one for half its mass, whose
   maximum differs, so the reported intensity is scaled by the wrong factor.
5. **`aborts_` is written from inside the parallel region.** `abort_`
   (`.cpp:1129-1140`) increments `aborts_[reason]` and, in debug mode, assigns
   `abort_reasons_[seed]`, both inside the `omp parallel for` of step 3.3 with
   no synchronisation. Two threads aborting at once corrupt the `std::map`; the
   C2 driver therefore records the library map single-threaded only.

### From the seed stage (B6)

Items 1 and 2 are executed; the others come from source review.

1. **`mass_trace:min_spectra = 1` silently finds nothing.** `min_spectra_ =
   floor(1 * 0.5) = 0`. Every trace score becomes 0/0 (`.cpp:340`), every peak a
   local maximum, and every overall score NaN. The run reports 0 seeds and
   0 features without a warning (executed: B6 driver `ffc1_min_spectra_1`).
   The parameter's minimum should be 2, or the value should be rejected.
2. **The overall score depends on the platform `powf`.**
   `std::pow(float, float)` at `.cpp:506` calls the C library `powf`, which is
   not correctly rounded on macOS: 99 of 30,840 executed scores are one binary32
   step off, so a score next to `seed:min_score` can flip a seed between
   platforms. Evaluating in `double` and rounding once gives the correctly
   rounded value.
3. **`write_debug` throws.** `writeFeatureDebugInfo_` reads
   `debug:pseudo_rt_shift` (`.cpp:2137`), but the declared parameter is
   `advanced:pseudo_rt_shift` (`.cpp:124`). The resulting `ElementNotFound`
   escapes the OpenMP region.
4. **Zero-width input ranges are undefined behaviour.** With one RT or one m/z,
   `intensity_rt_step_` or `intensity_mz_step_` is 0, and `intensityScore_`
   converts `floor(NaN)` to `UInt` (`.cpp:1837-1838`).
5. **`charge_low > charge_high` wraps.** `UInt charge_count = charge_high -
   charge_low + 1` (`.cpp:197`) wraps for `charge_low > charge_high + 1`, and
   the float arrays are then indexed past their end.
6. **`isotopeScore_` narrows its candidate search.** The inner loop starts at
   `best_end` (`.cpp:1759`), which a better fit found earlier in the outer loop
   has already raised. Combinations with fewer trailing isotopes are then never
   tried for later `b`, even when they fit better
   (`isotope_score_narrows_later_candidates_after_a_new_best_fit`).
7. **The peak-count check has a ranges message.** `getSize() == 0` reports
   "FeatureFinder needs updated ranges on input map" (`.cpp:1069-1071`), which
   describes neither the check nor the cause.
8. **The abundance override keeps a stray `(0, 1)` peak** (`.cpp:163-179`).
   Recorded by B2 in [ISOTOPE_SOURCE_PRECISION_SUPPORT](ISOTOPE_SOURCE_PRECISION_SUPPORT.md);
   C2 shows its effect on FFC_1: 0 seeds at 12C = 90 %.

## Ledger notes

- `FeatureFinderAlgorithmPicked.h` can move from `partial` to `complete` for the
  algorithm itself: every public and protected member of the header is ported,
  documented and tested, and `run` produces features end to end. Two members are
  deliberately not ported and are recorded as such in the API mapping:
  `writeFeatureDebugInfo_` and `abort_reasons_`, both reachable only through
  `write_debug`, which the port refuses because the source throws there. The
  `FeatureFinderDefs` struct in the same header is not ported and not used.
  `evidence_requires_review` remains the honest default until the integrator has
  read the mapping.
- Rust files: `src/analysis/feature_finder_picked/algorithm.rs`, `scoring.rs`,
  `seeds.rs`, `extension.rs`, `fitting.rs` and `resolution.rs`. Tests:
  `tests/feature_finder_picked_seeds.rs` and `tests/feature_finder_picked.rs`.
- The two behaviours B6 flagged are settled by the lead decision of 2026-09-15:
  `mass_trace:min_spectra = 1` follows the source (`CPP-271`), and a changed
  abundance computes the intended override (`CPP-247`). The second is a
  deliberate difference from the executed C++ and is the only one that changes
  which features are found.
- `docs/doc-coverage.json` needs `--write`: the six modules of this group are at
  100 %, and the floor was not re-recorded here because the file is the
  integrator's.
