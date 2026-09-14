# FeatureFinderAlgorithmPicked

This document covers the port of `FEATUREFINDER/FeatureFinderAlgorithmPicked.h`
and its implementation at core SDK `bc9cc12514c768385ce121d6ca4bb710fe1983c4`.
The source class is the centroided-peptide feature finder that
FeatureFinderCentroided runs.

The port is staged. Package B6-FFAP-SEEDS of the early TOPP bundle ports the
front half: the parameters, the input checks of `run`, the intensity, trace and
isotope-pattern scores, the isotope-pattern precalculation and seed selection
(steps 0 to 3.2 of `run_`). Seed extension, trace fitting, the quality checks
and overlap resolution (step 3.3 onward) belong to B7. Until then `run` returns
`Error::Unsupported` after selecting seeds.

| Rust file | Content |
| --- | --- |
| [`algorithm.rs`](../src/analysis/feature_finder_picked/algorithm.rs) | defaults, `Settings`, `validate_input`, `Limits`, `Options`, `run` |
| [`scoring.rs`](../src/analysis/feature_finder_picked/scoring.rs) | `position_score`, `nearest_from`, `IntensityThresholds`, `ScoreArrays`, `find_isotope`, `isotope_score` |
| [`seeds.rs`](../src/analysis/feature_finder_picked/seeds.rs) | `IsotopeWindows`, `SeedStage`, `ChargeSeeds`, `UserSeed`, `overall_score` |

The helper types come from
[`helper_structs.rs`](../src/analysis/feature_finder_picked/helper_structs.rs)
([support](FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md)). The binary32
isotope patterns come from `src/chemistry/isotopes.rs`
([support](ISOTOPE_SOURCE_PRECISION_SUPPORT.md)). Neither file is edited here.

Tests: [`tests/feature_finder_picked_seeds.rs`](../tests/feature_finder_picked_seeds.rs).
Manifest: [`tests/data/feature_finder_picked_provenance.json`](../tests/data/feature_finder_picked_provenance.json).
Nothing is feature-gated in the library. The integration test needs `mzml` and
`paramxml`; its user-seed case also needs `featurexml`.

**Module edges.** The package adds two cross-module edges:
- `analysis -> math`, for `pearson_correlation_coefficient`;
- `analysis -> param`, for `Param` and `DefaultParamHandler`.

Both are acyclic; `tools/check_module_cycles.py` accepts them.

## API mapping

Every member of the header is listed. "B7" marks members of the back half,
which that package ports.

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
| `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap&)` | `run`, `run_with_options` | returns `Error::Unsupported` after step 3.2 until B7 |
| `getDefaultParameters()` | `default_parameters()` | 29 entries, 7 section descriptions |
| `run_()` (protected) | `SeedStage::compute` for steps 0 to 3.2 | steps 3.3 and 4: B7 |
| `map_` | `SeedStage::experiment`, `into_experiment` | |
| `features_` | `RunOutput::features` | |
| `log_`, `debug_` | `Settings::write_debug` | debug output refused (native differences) |
| `aborts_`, `abort_reasons_`, `abort_()` | not ported | B7 |
| `seeds_` | `SeedStage::user_seeds` (`UserSeed`: m/z and RT) | the only fields the source reads |
| `pattern_tolerance_` ... `max_feature_intersection_`, `reported_mz_` | `Settings` fields of the same names; `reported_mz` as `ReportedMz` | |
| `intensity_rt_step_`, `intensity_mz_step_`, `intensity_thresholds_` | `IntensityThresholds::rt_step`, `mz_step`, `quantiles` | |
| `isotope_distributions_` | `IsotopeWindows::patterns` | |
| `updateMembers_()` | `Settings::from_parameters` | also reads the `run_` locals: charges, `fit:max_iterations`, abundances, seed thresholds, user-seed tolerances, `feature:min_score`, `feature:rt_shape` |
| `intersection_()` | not ported | B7 |
| `getIsotopeDistribution_(double)` | `IsotopeWindows::get` | |
| `findBestIsotopeFit_()` | not ported | B7 |
| `extendMassTraces_()`, `extendMassTrace_()` | not ported | B7 |
| `nearest_()` | `nearest_from` | also returns the steps walked |
| `findIsotope_()` | `find_isotope` | returns work units; checks indices |
| `positionScore_()` | `position_score` | |
| `isotopeScore_()` | `isotope_score` | |
| `intensityScore_(Size spectrum, Size peak)` | `IntensityThresholds::score(rt, mz, intensity)` | |
| `intensityScore_(Size rt_bin, Size mz_bin, double)` | `IntensityThresholds::bin_score` | `None` outside the grid |
| `chooseTraceFitter_(double&)` | `Settings::rt_shape` (`RtShape`) | the fitters: B4, B5, B7 |
| `cropFeature_()`, `checkFeatureQuality_()` | not ported | B7 |
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
`IsotopeWindows::mass_window_width`, and the `SeedStage` accessors.

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

## Native differences

| Source behaviour | This port | Reason and evidence |
| --- | --- | --- |
| `run` continues with seed extension, fitting and overlap resolution | returns `Error::Unsupported` after step 3.2 | B7 ports the back half. |
| scores are float data arrays appended to each spectrum, replacing its existing float arrays | `ScoreArrays`: one flat `f32` array per score, outside the spectra | Only the algorithm reads them, and only debug mode writes them. The input spectra keep their arrays, and no per-spectrum allocation is needed. |
| `spectrumRanges().byMSLevel(1)` needs `updateRanges()` from the caller | ranges are computed on demand from the validated spectra | No stale-range state exists, so the source's FAIMS "No ranges for this MS level" crash cannot occur. The source message "needs updated ranges" belongs to the peak-count check and is kept verbatim. |
| NaN or infinite RT, m/z or intensity values are sorted and binned with undefined results | `Error::InvalidValue` | The native readers never produce them. |
| `mass_trace:min_spectra = 1` gives `min_spectra_ = 0`. Every trace score becomes 0/0 = NaN and every peak a local maximum; no overall score reaches a threshold, no seed is found and the run returns an empty map | `Error::InvalidValue` | The plan's acceptance criterion requires an explicit error. The execution (B6 driver, `ffc1_min_spectra_1`) shows that the source is defined here: NaN trace scores, 0 seeds, 0 features, exit 0. The integrator decides whether the port follows the source instead. |
| a single retention time or a single m/z makes the bin width 0; `(UInt)floor(NaN)` is undefined | `Error::InvalidValue` | Undefined behaviour. |
| `charge_low > charge_high + 1` wraps the `UInt` charge count and indexes past the score arrays | `Error::InvalidValue` from `Settings::charge_count` | Undefined behaviour. `charge_low == charge_high + 1` gives zero charges, as in the source. |
| `write_debug = true` writes `debug/` into the working directory and then throws on the undeclared `debug:pseudo_rt_shift` | `Error::Unsupported` | A defect, and library code does not write files. |
| a changed `abundance_12C` or `abundance_14N` builds the override from a default `IsotopeDistribution` that already holds `(0, 1)`; the patterns grow (FFC_1 window 0: 27 normalised bins) and FFC_1 with 12C = 90 % finds 0 seeds (C2 `ffap_ffc1_abundance_12C_90`) | `Error::Unsupported` by default; `AbundanceOverride::Intended` uses the intended two-isotope distribution | The generator rejects the stray-peak distribution, so the defect cannot be reproduced. The intended distribution matches the executed `set()` construction (B2 probe `set_12C_90`). Which behaviour is right is a scientific decision. |
| the overall score is `std::pow(float, float)`, the platform `powf` | `libm::pow` in `f64`, rounded once to `f32` | Apple `powf` misrounds 99 of the 30,840 retained overall scores by one binary32 step, 12 of 3,084 on FFC_1. The port's value is the correctly rounded power for all of them (60-digit decimal check) and is the same on every machine. `libm::powf` differed on 226 of 3,084. No seed list changes. `tests/data/feature_finder_picked/overall_rounding.tsv` lists every difference. |
| `std::sort` of seeds with equal `f32` intensity is unspecified | stable: scan, then peak order | Deterministic. No retained configuration has a tie (the drivers check adjacent ties). |
| `std::sort` of a bin's intensities: `-0.0` and `+0.0` compare equal | `total_cmp`: `-0.0` first | Only a quantile's zero sign can differ, and only for inputs with both signed zeros; the tool path filters non-positive intensities. |
| progress and `Found N seeds for charge c.` go to `std::cout`; warnings to `OPENMS_LOG_WARN` | `SeedStage::log`, `RunOutput::log` | Library code never prints. B7 must interleave its `Found N feature candidates` lines per charge. |
| steps 3.1 to 3.3 run per charge | steps 3.1 and 3.2 run for every charge first | Step 3.3 reads only its own charge's arrays, so the arrays and seeds are identical. |
| user seeds are a copied `FeatureMap` sorted with `std::sort` | positions only, sorted stably; non-finite positions are `Error::InvalidValue` | NaN breaks the source sort. Equal m/z values give the same search result in any order. |
| unbounded work | `Limits`: spectra, peaks, charges, bins per dimension, windows, pattern values, score bytes, work units | Checked before the allocation or computation each bounds. The FFC_1 workload is several orders of magnitude below every default. |
| `Math::pearsonCorrelationCoefficient` returns an infinity when its denominator underflows to 0 from non-zero deviations | NaN, which counts as 0 | Inherited from `src/math/statistic_functions.rs`. Unreachable at isotope intensity scales. |
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
is compared bit for bit on Linux x86-64. The only exception is the documented
overall-score rounding, whose correctly rounded values the test asserts instead.

| Evidence | Configurations | Test |
| --- | --- | --- |
| C1 `FFC_write_ini`: FeatureFinderCentroided `-write_ini` (two runs, identical) | the whole algorithm section | `default_parameters_equal_the_executed_write_ini_section` |
| C2 `ffap_stages` library state (`omp1` and `omp4`, identical): effective members, bins, 21 quantiles, 15 windows, all score arrays of 3,084 peaks, printed seed counts | FFC_1 INI (tool loading), class-test INI, #9247 tight pattern and tight trace, FFC_1 with the retained output as user seeds | `ffc_1_stage_...`, `class_test_stage_...`, `tolerance_swap_stages_...`, `user_seed_stage_matches_the_executed_library` |
| B6 `seed_stage` library state (two runs, byte-identical), same content | default parameters (10 bins, charges 1 to 4, width 25, 107 windows); FFC_1 INI with 7 bins and charges 1 to 3 (21 windows) | `default_parameter_stage_...`, `seven_bin_three_charge_stage_matches_the_executed_library` |
| B6 `seed_stage`, `ffc1_min_spectra_1` | `mass_trace:min_spectra = 1`: `min_spectra_` 0, 0 seeds, 0 features, exit 0 | `min_spectra_one_is_refused` |

The FFC_1 score table also pins every loaded peak's m/z and intensity bits
against the C++ loader.

The oracle has two single-bin configurations (FFC_1 and the class test). The
B6 driver adds 7 and 10 bins, which exercise the four-cell interpolation, and up
to four charges.

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

- **FeatureFinderAlgorithmPicked_test.** The constructor and destructor
  sections are covered by `constructor_and_destructor`. The class test has no
  `getDefaultParameters` section of its own; the executed `-write_ini`
  comparison covers it. The `run()` and `[EXTRA] #9247` sections produce
  features and belong to B7; their seed stages are tier 1 above.
- **Source review.** Covered cases:
  - the input check order, including invalid parameters on empty input and on
    MS2 input;
  - unsorted input: sorted with the warning, and scored exactly like the sorted
    input;
  - the refusals: `write_debug`, charge wrap, zero-width ranges, non-finite
    user seeds, abundances;
  - the conversions in `settings_follow_update_members`;
  - restriction and type violations.
- **Hand-derived.**
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
- **Limits.** Each ceiling fails one below the FFC_1 requirement and passes at
  it (`limits_are_checked_before_the_work`).

### Verification commands

Run from the worktree root on the IBMI build node through
`~/.local/bin/openms-kim-gate.sh`; `target_verification` in the manifest records
the results.

```text
cargo fmt --all -- --check
cargo test --locked --all-features --test feature_finder_picked_seeds \
  --test feature_finder_picked_helper_structs --test isotopes_source_precision
cargo +1.85.0 test (same targets)
cargo test --locked --no-default-features --features mzml,paramxml,featurexml --test feature_finder_picked_seeds
cargo clippy --locked --all-features --all-targets -- -D warnings
cargo doc --locked --all-features --no-deps
python3 tools/check_core_sdk.py
python3 tools/check_module_cycles.py
python3 tools/check_doc_coverage.py
```

## C++ issue candidates

Candidates for the integrator's `OpenMS_CPP_ISSUES.md`. Items 1 and 2 are
executed; the others come from source review.

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

- `FeatureFinderAlgorithmPicked.h` is `partial`. The parameter surface, the
  input checks and the members of steps 0 to 3.2 are ported, documented and
  tested. Everything marked B7 is not ported.
- Rust files: `src/analysis/feature_finder_picked/algorithm.rs`, `scoring.rs`
  and `seeds.rs`. Tests: `tests/feature_finder_picked_seeds.rs`.
- The integrator decides on the flagged behaviours: `min_spectra = 1` and
  `AbundanceOverride`.
