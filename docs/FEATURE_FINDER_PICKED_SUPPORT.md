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
resolution (steps 3.3 and 4). `run` now produces features end to end. The
lanes port/ffap-semantics (non-finite and degenerate input, `FeatureFinderDefs`)
and port/ffap-instrumentation (the reusable instance, debug mode, progress)
completed the header, and their merge port/ffap-complete applied the lead's
wave-5 decisions D1 to D9 in one combined fix round: every source sort in the
Release build's order, the Release build's `powf`, the area iterator's
drift-time filter, the step-2.5 `length_error`, `ChargedIndexSet` equality
and glibc's `-nan`. A second combined fix round closed the findings of its
two verifiers: the debug side effects of a run that fails before the seed
loop, the tool's exit for the `length_error`, 64-bit integer parameters as the
source narrows them, the SSE operand order of the gnuplot formulas, the debug
files of a step-3.3.5 termination, and the introsort on comparators that are
not strict weak orderings. A third combined fix round applied lead decisions
D10 to D12: both trace fitters call the reference build's glibc `exp` and
`log`, ported, so every fit is the Release build's on every platform; step 2.5
follows the source past windows whose binary32 bins all underflow and past a
NaN `intensity_percentage_optional`; every `UInt` wrap that writes out of
bounds is refused at the wrap and the step-1 progress range wraps as executed;
the correlations are `Math::pearsonCorrelationCoefficient`'s, division by a
zero denominator included; unsorted input with mis-sized data arrays gives the
source's `Exception::Precondition`; and a seed-loop refusal where the executed
process dies records that termination with the seed's log lines. The
multi-thread race on `aborts_`, `abort_reasons_` and `log_` is the one
accepted exception to D1 (D11). A fourth combined fix round stored the
source's non-finite widths and meta values. A fifth recorded the terminations
outside the seed loop (the step-4 charge remainder, the stale abort seed, the
wrapped score-array count) and modelled the never-closed `log_` stream across
the runs of an instance, so that every termination states the length at which
the executed process leaves `debug/log.txt`; it also applied the lead's
decisions D13 (below, *Ledger notes*).

| Rust file | Content |
| --- | --- |
| [`algorithm.rs`](../src/analysis/feature_finder_picked/algorithm.rs) | defaults, `Settings`, `validate_input`, `Limits`, `Options` with `AbundanceOverride` and `DegenerateBinStep`, `PseudoRtShiftKey`, `RejectedParameters`, `run` and `run_with_options` (a fresh map), `feature_stage` (stage level) |
| [`scoring.rs`](../src/analysis/feature_finder_picked/scoring.rs) | `position_score`, `nearest_from`, `IntensityThresholds`, `ScoreArrays`, `find_isotope`, `isotope_score`; the crate-private `x86_64` emulation of the Release build's `UInt` conversion and NaN rules |
| [`defs.rs`](../src/analysis/feature_finder_picked/defs.rs) | `FeatureFinderDefs`: `IndexPair`, `IndexSet`, `ChargedIndexSet`, `Flag`, `NoSuccessor` |
| [`seeds.rs`](../src/analysis/feature_finder_picked/seeds.rs) | `IsotopeWindows`, `SeedStage`, `ChargeSeeds`, `UserSeed`, `overall_score` |
| [`extension.rs`](../src/analysis/feature_finder_picked/extension.rs) | `OverallScores`, `find_best_isotope_fit`, `extend_mass_traces`, `extend_mass_trace` |
| [`fitting.rs`](../src/analysis/feature_finder_picked/fitting.rs) | `FittedModel`, `crop_feature`, `check_feature_quality`, `build_feature`, the abort reasons |
| [`resolution.rs`](../src/analysis/feature_finder_picked/resolution.rs) | `intersection`, `resolve_overlaps`, `annotate_apex` |
| [`instance.rs`](../src/analysis/feature_finder_picked/instance.rs) | `FeatureFinderAlgorithmPicked`, the source object with its state across runs (parameters, seeds, aborts, abort reasons, the log stream, the isotope windows, progress), and its `run` into the caller's map |
| [`debug.rs`](../src/analysis/feature_finder_picked/debug.rs) | `DebugOutput`, `DebugLog`, `SeedMap`, `AbortReasons`, `seed_map`, `abort_map`, `debug_experiment`, `write_feature_debug_info`, `PseudoRtShift`, `HEAP_ADDRESS_END`, `ReportLine` |
| [`source_sort.rs`](../src/analysis/feature_finder_picked/source_sort.rs) | `source_sort_by`, `source_sort_reversed_by`, `source_sort_permutation` (the Linux x86_64 Release build's `std::sort`, libstdc++'s introsort) and `source_stable_sort_permutation` with `TemporaryBuffer` (its `std::stable_sort`), comparison by comparison, NaN keys included |
| [`glibc_powf.rs`](../src/analysis/feature_finder_picked/glibc_powf.rs) (crate-private) | `powf` as the reference build's GNU C Library 2.39 computes it (`__powf_fma`), ported from Arm optimized-routines (MIT; the notice is in the file), and `mul`, SSE's `mulss` NaN rule |
| [`glibc_libm.rs`](../src/analysis/feature_finder_picked/glibc_libm.rs) (crate-private) | `exp` and `log` as the reference build's GNU C Library 2.39 computes them (`__ieee754_exp_fma`, `__ieee754_log_fma`), ported from Arm optimized-routines (MIT; the notice is in the file); `atan` (the host's with glibc, the `libm` crate's elsewhere, lead decision D10's fallback) and `sqrt` with SSE's NaN bits |

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
Of those files the combined fix rounds changed only what this algorithm
reaches: `trace_fitter::stream_number` prints glibc's `-nan`; both fitters'
`gnuplot_formula` computes its sums and products in the Release build's SSE
operand order; both fitters call `glibc_libm` instead of the platform's or the
`libm` crate's `exp`, `log`, `atan` and `sqrt` (round 3); `helper_structs.rs`
names its NaN-merge refusal text; and `isotopes.rs` gains the crate-private
`CoarseIsotopePatternGenerator::estimate_from_peptide_weight_source`, which
returns `SourceSingleEstimate::AllUnderflowed` where every retained binary32
bin underflows instead of the error its public callers keep.

Tests: [`tests/feature_finder_picked_seeds.rs`](../tests/feature_finder_picked_seeds.rs)
(steps 0 to 3.2, the degenerate intensity bins and `FeatureFinderDefs`),
[`tests/feature_finder_picked.rs`](../tests/feature_finder_picked.rs)
(steps 3.3 and 4, the intended abundance override, the non-finite, sort,
drift-time and boundary captures, the charge-count wraps and every source sort
against the executed library), the unit tests of `glibc_powf.rs`,
`glibc_libm.rs` and `source_sort.rs` (`--lib`) and, through the tool,
[`tests/topp_feature_finder_centroided.rs`](../tests/topp_feature_finder_centroided.rs);
[`tests/feature_finder_picked_instrumentation.rs`](../tests/feature_finder_picked_instrumentation.rs)
(the instance, debug mode, progress, a caller's map, 64-bit integer
parameters, a failed step 2.5, the gnuplot formulas on non-finite
operands, the seed-loop crashes and the step-1 progress wrap; manifest
[`tests/data/feature_finder_picked_instrumentation_provenance.json`](../tests/data/feature_finder_picked_instrumentation_provenance.json),
which needs `featurexml` too).
Manifest: [`tests/data/feature_finder_picked_provenance.json`](../tests/data/feature_finder_picked_provenance.json).
Nothing is feature-gated in the library. Both integration tests need `mzml` and
`paramxml`; the user-seed cases also need `featurexml`.

**Module edges.** The packages add three cross-module edges, all of which
the recorded graph already held, so `tools/check_module_cycles.py` reports no
new edge:
- `analysis -> math`, for `pearson_correlation_coefficient` (whose range
  error `scoring::source_pearson` returns);
- `analysis -> param`, for `Param` and `DefaultParamHandler`;
- `analysis -> concept`, for `PROTON_MASS_U`, `UserParam::NUM_OF_DATAPOINTS`
  and `parallel::{Threads, map_collect}`.

`analysis -> kernel` and `analysis -> metadata` were already there; B7, the two
lanes and the combined fix round add no further edge (`source_sort.rs` and
`scoring.rs` use each other inside the module).

## API mapping

Every member of the header is listed.

### `FeatureFinderAlgorithmPicked`

| Source member | Rust | Notes |
| --- | --- | --- |
| base `DefaultParamHandler` | `instance::FeatureFinderAlgorithmPicked`: `name`, `set_name`, `subsections`, `handler_equal` (`operator==`), `defaults`, `default_parameters`, `parameters`, `set_parameters`, `set_parameters_logged`, over the crate's `DefaultParamHandler`; `Settings::from_parameters` for the stateless `run`. The static `writeParametersToMetaValues` is `DefaultParamHandler::write_parameters_to_meta_values`; the protected `check_defaults_` and `warn_empty_defaults_` are not settable, as the source object never changes them | `setName` renames the handler in the unknown-parameter warnings, `getSubsections` is empty, and `operator==` compares only the handler part, with the refused set the source keeps in `param_` (source-reviewed, `the_handler_base_renames_compares_and_has_no_subsections`). The checks are `Param::checkDefaults`'s, with its messages (`algorithm::check_parameters`): an integer is narrowed to its low 32 bits (`int tmp = value`) before its restriction is checked, so `intensity:bins = 2^32 + 10` is accepted and `2^32` refused as `'0'` (executed: `param_narrowing.tsv`, 21 cases, `integer_parameters_beyond_the_int_range_follow_the_release_build`). Executed: `getParameters` after construction, after a run and after an empty-input run; the unknown-parameter warnings in `checkDefaults` order, and on a refused set only those before the refused entry; after a refused set `parameters()` shows that set merged with the defaults while the settings stay, as the source assigns before it checks (`algorithm::RejectedParameters::Shown`, the default; `Discarded` keeps the accepted set) |
| base `ProgressLogger` | `instance::FeatureFinderAlgorithmPicked::set_log_type`, `log_type`, `set_progress_logger`, `progress_logger_mut`; the base's public `startProgress`, `setProgress`, `nextProgress` and `endProgress` are the methods of the logger `progress_logger_mut` returns, and are absent under type `NONE` (where the source's calls do nothing), as `progress_logger_mut` then returns `None` | the 20 call sites of `.cpp:241-992`: every start, set and end call with its label, range and value, equal to the Release build's complete event sequence (executed with a counting clock, see *Debug mode*); no logger (type `NONE`, the default) costs nothing. An inverted range (steps 2 and 3.2 on fewer than `2 * min_spectra` scans) is passed unchanged, as in the Release build |
| `MapType`, `SpectrumType`, `FloatDataArrays` | `MSExperiment`, `MSSpectrum`, `ScoreArrays` | see *Score arrays* below |
| `PeakType`, `Seed`, `MassTrace`, `MassTraces`, `TheoreticalIsotopePattern`, `IsotopePattern` (protected) | `Peak1D` and the `helper_structs` types | |
| `FeatureFinderAlgorithmPicked()` | `instance::FeatureFinderAlgorithmPicked::new`, `with_options`; `default_parameters`, `HANDLER_NAME` | |
| `setSeeds(const FeatureMap&)` | `instance::FeatureFinderAlgorithmPicked::set_seeds`, `seeds`; `run` replaces them, as the source's `run` does, and sorts them in place (`seeds_.sortByMZ()`, the Release build's introsort); the `seeds` argument of `SeedStage::run` | every run replaces them with the caller's map and sorts that, so `seeds()` returns the last run's seeds in sorted order (unsorted when `run_` threw before `:190`, executed: `run_max_iterations_negative`) |
| `setData_(MSExperiment&&, FeatureMap&)` (private) | `instance::FeatureFinderAlgorithmPicked::run` consumes the experiment and extends the caller's `&mut FeatureMap` | |
| `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap&)` | `instance::FeatureFinderAlgorithmPicked::run`, which extends the caller's map as the source does (see *Reusing an instance*); `algorithm::run` and `run_with_options` for a fresh map | the stateless functions return `RunOutput` with the features, the log, the abort counts and the debug output |
| `getDefaultParameters()` | `default_parameters()` | 29 entries, 7 section descriptions |
| `run_()` (protected) | `instance` `run_core`: `SeedStage::prepare` (steps 0 to 2.5), then per charge `select_next_charge`, `extend_charge` and `settle_charge`, then step 4; `SeedStage::compute` and `feature_stage` at stage level | |
| `map_` | `SeedStage::experiment`, `into_experiment` | |
| `features_` | the caller's `&mut FeatureMap` of `instance::FeatureFinderAlgorithmPicked::run`; `RunOutput::features` of the stateless `run` | |
| `log_`, `debug_` | `Settings::write_debug`; `debug::DebugOutput` (`log`, `log_opened`) returned by `run` | the log text in source order and formatting, with the stream's state across runs, also after a run that fails in steps 1 to 2.5 (see *Debug mode*); the tool writes the file |
| `aborts_`, `abort_()` | `instance::FeatureFinderAlgorithmPicked::aborts` (`u32`, accumulating over runs); `RunOutput::aborts` | counted serially in seed order: the source's single-thread counts at any thread count; the source races them (candidate 5) |
| `abort_reasons_` | `instance::FeatureFinderAlgorithmPicked::abort_reasons`, `debug::AbortReasons` (keyed by intensity, as `Seed::operator<`), `debug::abort_map` | never cleared, as in the source; the abort map reads the stored indices in the current input (see *Reusing an instance*) |
| `seeds_` | `instance::FeatureFinderAlgorithmPicked::seeds`; `SeedStage::user_seeds` (`UserSeed`: m/z and RT, the only fields the source reads) | sorted by `seeds::sort_user_seeds` (crate-private) |
| `pattern_tolerance_` ... `max_feature_intersection_`, `reported_mz_` | `Settings` fields of the same names; `reported_mz` as `ReportedMz` | |
| `intensity_rt_step_`, `intensity_mz_step_`, `intensity_thresholds_` | `IntensityThresholds::rt_step`, `mz_step`, `quantiles` | |
| `isotope_distributions_` | `IsotopeWindows::patterns`; `IsotopeWindows::precalculate_onto` and `instance::FeatureFinderAlgorithmPicked::isotope_windows` keep them across runs; `IsotopeWindows::source_count`, `SOURCE_MAX_WINDOWS`, `LENGTH_ERROR_WHAT` | extended, never cleared, as in the source (see *Reusing an instance*); `resize` above `vector::max_size()` is its `std::length_error` text |
| `updateMembers_()` | `Settings::from_parameters`; in the instance `Settings::update_members` (crate-private), in the source's order | the conversions are the Release build's: `operator unsigned int` for `max_missing` and `bins` (a negative value throws `ConversionError` half way, the earlier members keep their new values), and `min_spectra_` as the `cvttsd2si` of `floor(value * 0.5)`, low 32 bits (`libOpenMS.so` `0x18da803`-`0x18da866`). `Settings` also holds the `run_` locals (`read_run_values`, crate-private): charges (`(Int)`, low 32 bits), `fit:max_iterations` (`operator unsigned int`, whose `ConversionError` for a negative value the run reports where `run_` throws it, before the user seeds are sorted), abundances, seed thresholds, user-seed tolerances, `feature:min_score`, `feature:rt_shape` |
| `intersection_()` | `resolution::intersection` | |
| `getIsotopeDistribution_(double)` | `IsotopeWindows::get` | |
| `findBestIsotopeFit_()` | `extension::find_best_isotope_fit` | returns the score and the pattern instead of an out-parameter |
| `extendMassTraces_()`, `extendMassTrace_()` | `extension::extend_mass_traces`, `extension::extend_mass_trace` | the overall-score array is `extension::OverallScores` |
| `nearest_()` | `nearest_from` | also returns the steps walked |
| `findIsotope_()` | `find_isotope` | returns work units; checks indices |
| `positionScore_()` | `position_score` | |
| `isotopeScore_()` | `isotope_score` | |
| `intensityScore_(Size spectrum, Size peak)` | `IntensityThresholds::score(rt, mz, intensity)` | total: the undefined `UInt` conversion and the NaN bits follow the Linux x86_64 Release build (*Degenerate intensity bins*) |
| `intensityScore_(Size rt_bin, Size mz_bin, double)` | `IntensityThresholds::bin_score` | `None` outside the grid; a NaN result carries the Release build's bits |
| `chooseTraceFitter_(double&)` | `fitting::FittedModel::new` from `Settings::rt_shape` | the enum replaces the pointer plus the `tau != 0` flag and the `dynamic_pointer_cast`. Its two `OPENMS_LOG_DEBUG` lines (`use asymmetric rt peak shape`, `use symmetric rt peak shape`, `.cpp:1902`, `:1908`) are not ported: they go to the thread-local debug stream (`getThreadLocalLogDebug`), which the crate has no sink for, and the executed FeatureFinderCentroided prints none of them at `-debug 1` or `-debug 5`, on stdout or in `TOPP.log` (`../oracle/ffc-instrumentation-v4`, `tool_debug1`, `tool_debug5`, `tool_debug1_egh`), so no tool output depends on them; a library caller reads the choice from `Settings::rt_shape` |
| `cropFeature_()` | `fitting::crop_feature` | returns the cropped traces instead of an out-parameter |
| `checkFeatureQuality_()` | `fitting::check_feature_quality` | returns `QualityOutcome`: `Accepted(FeatureQuality)` or `Rejected(reason)` |
| step 3.3 (`.cpp:576-856`), the `omp parallel for` and the containment pass | `instance::extend_charge` and `settle_charge`, with `fitting::build_feature` for step 3.3.5 | `concept::parallel::map_collect` over the seed indices |
| step 4 (`.cpp:859-1016`) | `source_sort::source_sort_by` (the Release build's introsort order), `resolution::resolve_overlaps`, `retain`, `source_sort_by` with the arguments swapped, `resolution::annotate_apex` | |
| `writeFeatureDebugInfo_()` | `debug::write_feature_debug_info`, `debug::FeatureDebugFiles`, `debug::PseudoRtShift`; `algorithm::PseudoRtShiftKey` | the `.dta`, `_cropped.dta` and `.plot` texts, byte-identical to the Release build; the undeclared `debug:pseudo_rt_shift` is read as the source reads it, a string or list value included (see *Debug mode*) |
| `operator=`, copy constructor (private, not implemented) | not applicable | the stage is an owned value |
| step 1, second half (`.cpp:280-287`) | `fill_intensity_scores` (crate-private) | |
| step 2 (`.cpp:291-348`) | `fill_trace_scores` (crate-private) | |
| step 3.1 (`.cpp:447-488`) | `fill_pattern_scores` (private) | |
| step 3.2 (`.cpp:489-574`) | `select_seeds` (private), `overall_score`, `SeedStage::charges` | `overall_score` computes the reference build's `powf` (`glibc_powf`) |

Native additions: `Limits`, `Options`, `AbundanceOverride`, `DegenerateBinStep`, `RunOutput`,
`validate_input`, `UNSORTED_WARNING`, `BASE_MAX_ISOTOPES`, the `source_sort`
functions and `TemporaryBuffer`, `IsotopeWindows::source_count`,
`SOURCE_MAX_WINDOWS`, `LENGTH_ERROR_WHAT`,
`OVERRIDE_EXTRA_ISOTOPES`, `Settings::charge_count`, `Settings::max_isotopes`,
`QUANTILE_COUNT`, `ScoreArrays`, `ChargeSeeds`, `UserSeed`,
`IsotopeWindows::mass_window_width`, the `SeedStage` accessors, and from the
back half `OverallScores`, `FittedModel`, `FeatureQuality`, `QualityOutcome`,
`FeatureInput`, the seven `ABORT_*` reason constants, `SPECTRUM_INDEX`,
`SPECTRUM_NATIVE_ID` and `invalid_apex_warning`.

### `FeatureFinderDefs` (same header)

| Source member (`.h:24-55`) | Rust (`defs.rs`) | Notes |
| --- | --- | --- |
| `IndexPair` (`IsotopeCluster::IndexPair`, `std::pair<Size, Size>`) | `IndexPair = (usize, usize)` | the source comment says two `UInt`s; they are `Size` |
| `IndexSet` (`std::set<IndexPair>`) | `IndexSet = BTreeSet<IndexPair>` | the same lexicographic order |
| `ChargedIndexSet` (derived from `IndexSet`, `Int charge`, constructor sets 0) | `ChargedIndexSet { indices, charge }`, `Default`, `Deref`/`DerefMut` to the set | inheritance becomes a field plus `Deref` |
| `a == b`, `a < b`, ... on `ChargedIndexSet` (no operator of its own: `std::set`'s through the base class) | `PartialEq`, `Eq`, `PartialOrd`, `Ord` on `indices` only | the charge takes no part, as in the source (executed: `defs_eq_probe`, seven pairs, all six operators) |
| `enum Flag { UNUSED, USED }` | `#[repr(i32)] enum Flag { Unused = 0, Used = 1 }` | the Release build's enum is 4 bytes |
| `NoSuccessor(file, line, function, index)` | `NoSuccessor::new(index)` with `#[track_caller]` | the caller's file and line are recorded (`file()`, `line()`); Rust has no function name to record |
| `NoSuccessor::index_` (protected) | `NoSuccessor::index()` | |
| name `"NoSuccessor"`, `what()` | `NoSuccessor::NAME`, `name()`, `message()`, `Display` | `there is no successor/predecessor for the given Index: <scan>/<peak>` |
| `GlobalExceptionHandler::setMessage(what())` in the constructor | not ported | the port keeps no process-wide exception state |
| `~NoSuccessor()` | nothing | the default destructor |
| (exception type) | `std::error::Error`; `From<NoSuccessor> for Error` gives `Error::InvalidValue("NoSuccessor: <message>")` | a module-local error; the crate has no dedicated variant |

The algorithm uses none of them. The struct duplicates
`FEATUREFINDER/FeatureFinderDefs.h`, which nothing includes; a translation unit
that includes both headers does not compile (candidate 6 of the section below;
[helper structs support](FEATURE_FINDER_PICKED_HELPER_STRUCTS_SUPPORT.md)).
Evidence: `feature_finder_defs_match_the_executed_probe` against the executed
probe `defs_probe` (enumerator values, enum size, a default
`ChargedIndexSet`, set order, and name and message for three index pairs up to
`SIZE_MAX`), and `charged_index_set_comparisons_match_the_executed_probe`
against `defs_eq_probe` (`../oracle/ffap-complete-fix1`, two repetitions,
identical).

## Debug mode

`write_debug = true` makes the source write into the working directory
(`FeatureFinderAlgorithmPicked.cpp:226-232`, `550-571`, `714`, `1028-1052`).
The library port returns the same content as data in `debug::DebugOutput`
(`RunOutput::debug`, `FeatureFinderAlgorithmPicked::debug_output`), and
`FeatureFinderCentroided` creates `debug/` and `debug/features/` and writes the
files under the source's names:

| Source output | `DebugOutput` | Evidence (Linux x86_64 Release, `OMP_NUM_THREADS=1`) |
| --- | --- | --- |
| `debug/log.txt` (`log_`, 58 write statements) | `log` (`DebugLog`), `log_opened` | byte-identical for the tool cases a1, a2, a3 and the driver's two-run object; `double` values printed as `operator<<` prints them, glibc's `nan`/`-nan` included; a run that fails in step 2.5 leaves the first line (`../oracle/ffap-complete-fix2`, `lenerr_*` and `tool_*`, see below) |
| `debug/seeds_<charge>.featureXML`, per charge, also for a charge without seeds | `seed_maps` (`SeedMap`) | D6 (decoded, ids excluded): a1, a2, a3, a4, b1, stale scaled, debug_twice |
| `debug/features/<plot_nr>.dta`, `_cropped.dta`, `.plot` (`writeFeatureDebugInfo_`) | `feature_files` (`FeatureDebugFiles`, `debug::write_feature_debug_info`) | byte-identical: the 75 files and the log of each of the driver cases declared-shift500, -shift123, -int250, -egh and -prefilled, of a string and a string-list shift, and of a string shift with a scan at RT `5e-275` and `1e-289` (`debug_digests.tsv`); `double` values in the `.plot` formulas as `operator<<` prints them, glibc's `-nan` included, and `k * shift + rt` with x86_64's NaN rules (`inf * 0` is the negative default NaN on every host): a `+inf`, `-inf`, negative-NaN and positive-NaN shift (`shift_nonfinite_digests.tsv`, `../oracle/ffap-complete-fix1/node/run_shift.sh`, two runs each); the formulas' own sums and products follow the SSE operand order of `getGnuplotFormula` (`rt_shift` and `theoretical_int` are the destinations; EGH's `2 * sigma * sigma` is `(sigma + sigma) * sigma`), so a NaN the formula's own arithmetic creates, or passes on from its operands, prints with the executed sign on every host (`gnuplot_formula_nonfinite.tsv`, 648 executed formulas, `gnuplot_formulas_print_the_executed_nan_signs`). A NaN operand the fit produced keeps the sign the fit gave it; see *The sign of a NaN inside a fit* below |
| `debug/abort_reasons.featureXML` | `abort_reasons` (`debug::abort_map`) | D6: a1, a2, declared-*, debug_twice, stale scaled; the feature ids `0, 1, ...` exactly |
| `debug/input.mzML`: the input with the score arrays, without the overall score | `input` (`debug::debug_experiment`) | a1, a2, a3, debug_twice: every float array bit for bit, NaN bits and overall scores included. `mzml::write_source_float_arrays` writes the non-finite values the source writes |
| the process terminates: in the seed loop in `writeFeatureDebugInfo_`, after it at step 3.3.5 (`.cpp:790`), at the out-of-bounds read of an empty best pattern in `extendMassTraces_`, or never returns from the NaN profile merge (`CPP-242`); before it at a wrapped score-array count (`.cpp:196-221`); after it at the step-4 charge remainder (`.cpp:936`, `:945`) or a stale abort seed (`.cpp:1037-1039`) | `termination` (`DebugTermination`: its `TerminationPoint` `Seed`, `ScoreArrays`, `OverlapResolution` or `AbortMap`, its `TerminationKind` `Exception`, `OutOfBounds`, `ArithmeticTrap` or `NeverReturns`, and `log_file_bytes`, the length the executed process leaves `debug/log.txt` at), with `Error::Unsupported` or `Error::InvalidValue`; `FeatureFinderAlgorithmPicked::termination` records the same for a run without `write_debug` | a4 and b1: charge, exception and message; `log_file_bytes` is the length of the executed file after the SIGABRT. The terminations outside the seed loop and across runs: *The file at termination* below. At step 3.3.5 the seed's log lines and feature files are kept before the termination is recorded (source review; no executed input reaches it; `a_step_3_3_5_termination_keeps_the_seed_debug_output`). At an empty best pattern the refused seed's lines are appended first, and the executed `debug/log.txt` after the SIGSEGV is exactly the flushed prefix, with every executed feature file byte for byte (`neg_oob1`, `neg_oob_seed035`, `neg_none_avg0`: 53, 51 and 26 plots, `crash_digests.tsv`, `a_seed_loop_crash_keeps_what_the_executed_process_had_written`) |

**The undeclared key.** `writeFeatureDebugInfo_` reads
`param_.getValue("debug:pseudo_rt_shift")` (`.cpp:2137`), a key the defaults
do not declare; the declared one is `advanced:pseudo_rt_shift`. Without the key
`Param::getValue` throws `ElementNotFound`, inside the `omp critical` section of
the parallel seed loop, so the exception leaves the OpenMP region and
`std::terminate` ends the process: the executed FeatureFinderCentroided prints
OpenMS's fatal-exception block and dies of SIGABRT (shell status 134, cases a4
and b1, all repetitions). Every debug run in which a seed reaches the fit ends
this way. `Options::pseudo_rt_shift` chooses the port's behaviour:

- `PseudoRtShiftKey::Source` (the default, and the tool's): read
  `debug:pseudo_rt_shift` as the source does. An integer or float value is
  used; a missing key (`ElementNotFound`) or an empty value
  (`ConversionError`) stops the run at the first seed that reaches the fit
  with `Error::Unsupported`, after everything the source did before that
  point. A string or list value does not throw: `ParamValue::operator
  double()` returns the union member `dou_` for every type but `EMPTY` and
  `INT` (one `movsd 0x8(%rdi),%xmm0` in the Release `libOpenMS.so`), so the
  shift is the bit pattern of a heap pointer (`pun_values.txt`: nine
  values in three processes, all different). A heap pointer of a Linux
  x86_64 process lies below `DEFAULT_MAP_WINDOW = 2^47 - 4096`
  (`debug::HEAP_ADDRESS_END`; the reference host has 48-bit virtual
  addresses and no `la57`, and its kernel headers are recorded in
  `shift_band/address_space.txt`), so the shift is a subnormal number below
  `2^-1027`. Trace `k` writes `k * shift + rt` in the `.dta` files and
  `k * shift + centre` with six digits in the `.plot` file. Where the largest
  possible shift, `(2^47 - 4096) * 2^-1074 * k`, leaves that text as it is,
  every address does (the arithmetic and the rounding are monotone), and the
  files are those of shift 0: the Release build wrote the same 75 files and
  log in three processes each with a string and a list value, and with the
  first fitted seed's scan moved to RT `5e-275` and to `1e-289`
  (`ffap_shift_band_driver.cpp`), with a different address in every process,
  and the port writes those bytes (`debug::PseudoRtShift::HeapAddress`).
  Where that shift changes the text, the written number depends on the
  address: the scan moved to RT `0`, `1e-295` and `1e-300` wrote a different
  `0.dta` in each of three processes. The port refuses exactly there, at the
  first such value in the order the source writes them (`Error::Unsupported`
  at that seed's files), after the seed map and the log up to that point,
  which match the executed ones. For trace 1 the boundary lies at a positive
  retention time of `2^-974` (about `6.3e-294`; `2^-973` for trace 2); a
  fitted centre is refused below about `1.4e-303` (trace 1) whatever its
  digits, and up to about `2^-974 * k` (about `6.3e-294` for trace 1) when its
  six-digit text lies within `k * (2^47 - 4096) * 2^-1074` of a rounding
  boundary. A bound from the whole 56-bit address space would refuse `1e-289`,
  which the Release build writes reproducibly.
- `PseudoRtShiftKey::Declared`: read `advanced:pseudo_rt_shift` and write the
  feature files for every seed that reaches the fit, which is what the source
  evidently intends. The driver cases declared-* (with `debug:pseudo_rt_shift`
  set, so the Release build completes) give the same bytes under both
  policies.

A fit that throws `UnableToFit` (`TraceFitter.cpp:111`, `:129`) inside the
same parallel loop would end the source process at `.cpp:670`, before the
debug write at `:714`, but no input reaches either throw (see the
`Exception::UnableToFit` row of *Native differences*). The refused NaN profile
merge, where the source never returns, records a `NeverReturns` termination
after the seed's lines, like the out-of-bounds refusal above, with the plot
number the hanging seed received. Executed in fix round 4
(`../oracle/ffap-complete-fix4`, `fix4_vfi`, the round-3 instrumentation
verifier's driver; killed after 30 s, twice each, identical): with scan 50's
retention time NaN the Release process had flushed 1,016,234 bytes of
`debug/log.txt` and written the files of plots 0 to 2, Gaussian and EGH alike;
with scan 20's, 1,065,423 bytes and plots 0 to 9. The port's flushed prefix,
termination (plot numbers 3 and 10) and feature files equal them
(`a_seed_loop_that_never_returns_keeps_what_the_executed_process_had_written`,
`never_returns_digests.tsv`). An error from one
of the port's own ceilings ends the run with that error at the first such seed
in seed order and records no termination: the executed process would go on,
so the debug output holds the seeds before that one.

**A run that fails before the seed loop.** The source opens `log_` and
creates `debug/features` after the score arrays and before step 1
(`.cpp:226-232`) and writes `Precalculating intensity thresholds ...` at once.
A failure in steps 1 to 2.5 leaves both: at a maximum m/z of `1e19` step 2.5
throws `std::length_error`, the executed library's `debug/log.txt` is empty
while the object lives and holds those 40 bytes after it is destroyed, and
`debug/features` exists (`lenerr_single`); at `2e18` the allocation throws
`std::bad_alloc` with the same side effects. The port keeps them in
`DebugOutput` (the port's own window ceiling refuses the `2e18` count, after the
same point), and the stream stays open, so the object's next debug run writes
its seed map, abort map and input, byte-identical to a fresh object's in the
executed build, but no log (`lenerr_reuse`;
`a_debug_run_that_fails_in_step_two_point_five_keeps_the_executed_debug_output`).
A refusal before that point (the preflight ceilings, the charge wrap, the
score-array ceiling) leaves the stream closed, as the source has not reached
its `open`.

**The file at termination.** `log_` is an `std::ofstream` whose 8,191-byte
buffer is lost when the process aborts or crashes. `DebugLog` models
libstdc++'s `basic_filebuf` (`fstream.tcc`: a block write when the text does
not fit in the free space, one `sputc` per `char` insertion) and records how
many bytes reached the file. The model predicts the executed files of a4, b1,
of the two-run object while it was still alive and of the three SIGSEGV runs.
Because `log_` is an instance member that the first debug run opens and nothing
closes before the object is destroyed, the file a terminated process leaves is
the flushed prefix of whichever debug run of the object opened the stream:
this run's, or an earlier one's, whose complete log a caller wrote when that
run returned. Every `DebugTermination` therefore carries `log_file_bytes`, the
flushed length of that run (`None` when no debug run of the object opened the
stream), and `FeatureFinderAlgorithmPicked::debug_log_file` shows the stream
(`DebugLogFile`: written, flushed, failed). A caller writes `log.txt` in full
when a run that opened it returns, and cuts it to `log_file_bytes` at any
later termination; the tool, whose object runs once, writes that prefix.
Executed in fix round 5 (`../oracle/ffap-complete-fix5`, `fix5_driver` and the
unchanged `ffap_instr_driver`, every case twice, identical but for the abort
map's random unique id; `termination_digests.tsv.gz`):

- **Step 4.** A caller's charge-0 feature below the found charge-2 feature at
  m/z 652.766 makes step 4 compute `2 % 0`: SIGFPE (status 136), with and
  without `write_debug`. The debug run left `debug/log.txt` at 1,163,782 bytes
  (the `Intersection` line of the pair still in the buffer),
  `seeds_2.featureXML` and the files of 25 plots, and no abort map or input.
  The port refuses at that pair and records an `ArithmeticTrap` at
  `OverlapResolution`; its flushed prefix, seed map and 75 feature files are
  the executed ones. The same feature above the charge-2 one (`0 % 2 == 0`)
  and no extra feature return, with the complete log
  (`process_ending_refusals_outside_the_seed_loop_record_their_termination`).
- **The abort map.** A reused object whose abort seeds lie outside a
  four-scan input reads them after the second run's four seed maps: SIGSEGV
  (status 139), and `debug/log.txt` is the first run's flushed prefix,
  1,114,578 of its 1,118,230 bytes, as the second run's `open` failed. The
  port records an `OutOfBounds` termination at `AbortMap { entry: 0 }` with
  that length. On the same scans with doubled intensities the run returns
  and the destroyed object leaves the complete 1,118,230 bytes.
- **The score arrays.** `charge_low` 4 with `charge_high` 2 wraps the array
  count to one array, which the source writes past before it creates
  `debug/`: SIGSEGV and no file. The port records an `OutOfBounds` termination
  at `ScoreArrays` for the counts whose allocation it takes to succeed: the
  wraps to one array (`2^31 - 1`, `-1`) and, for `-4` and below, a wrapped
  count within its own ceiling `3 + 2 * Limits::max_charges`. Where the
  wrapped count is near `2^32` the executed allocation fails first and
  `std::bad_alloc` reaches the caller (7/2, below), which is no termination.
- **A later run of a reused object.** One object runs FeatureFinderCentroided_1
  with `write_debug` (the driver measured `debug/log.txt` at the flushed
  length while the object lived), then a second run with or without
  `write_debug`. Wherever the second run ends the process, the file is the
  first run's flushed prefix: 1,163,782 of 1,165,129 bytes after an empty best
  pattern in the seed loop (`vfi2_driver`'s `avg0` section, SIGSEGV) and after
  the wrapped score arrays 4/2 and `INT_MAX`/498 (SIGSEGV); 139,387 of 141,951
  bytes after a first run with a NaN `intensity_percentage_optional` and a
  charge-0 caller feature (SIGFPE). A plain first run leaves isotope windows
  that change the second run's features, so the charge-0 feature overlaps
  nothing and the run returns; the NaN cutoff empties the windows, and the
  second run computes a fresh object's. Where the second run returns (the
  same run again, or 7/2 whose `std::bad_alloc` the driver catches), the
  destroyed object leaves the complete first-run log. The port records every
  one of these terminations with the first run's length, none for the
  others, and the file protocol reproduces every executed file, feature files
  of both runs included
  (`a_reused_instance_leaves_the_executed_log_at_every_later_termination`).

**The sign of a NaN inside a fit.** The start values, the model, its bounds,
area and formulas, the cropping and quality scores and their correlations
compute every NaN they create with x86_64's rules on every host (the default
NaN `0xfff8000000000000` of an invalid operation, and the first NaN operand
otherwise), and the ported `exp` and `log` do too. The Levenberg-Marquardt
iterations themselves (`src/math/fitters/levenberg_marquardt.rs`) use the
host's arithmetic: a NaN they create from finite values (an overflowing step
followed by `inf - inf` or `0 * inf`) carries the host's sign, negative on
x86_64 and positive on arm64. No executed configuration reaches such a NaN;
if one did, the fitted value, and the `.plot` text or meta value that prints
it, would show `nan` on an arm64 host where the Release build shows `-nan`
(platform note, lead decision D8).

**Order.** The source writes the log from inside its parallel loop, under a
critical section, in schedule order, while `abort_` writes `aborts_`,
`abort_reasons_` and `log_` without one. The port collects every seed's lines
and aborts and appends them in seed order: the source's single-thread output at
every thread count (`debug_output_is_identical_at_every_thread_count`). The
executed tool at four threads wrote a different log in each of three
repetitions (cases c1, c2); that output is undefined and not compared. Lead
decision D11 accepts this as the one documented exception to D1: the race has
no reproducible result, and the port gives the single-thread result instead of
refusing, because the determinism contract requires parallel output to equal
serial output and a refusal would block every parallel run (`aborts_` is
written in every run).

**Progress.** The 20 `startProgress`/`setProgress`/`endProgress` calls
(`.cpp:241-992`) go to an optional `ProgressLogger` in source order. Two
executed comparisons pin them. The transcript with `setLogType(CMD)` matches
with the timing text masked (`the_progress_transcript_matches_the_release_build`);
it shows the labels and their order, but `setProgress` forwards a value only
when `time(nullptr)` has changed, so the values it prints depend on the wall
clock. The second driver (`ffap_progress_driver.cpp`) therefore defines
`time()` itself as a counter (exported with `-rdynamic`, so `libOpenMS.so`
binds to it) and installs a recording `ProgressLoggerImpl` with `setLogger`:
every call is forwarded and printed with its raw arguments, and the output is
identical in both repetitions. With a counting clock the port's logger
forwards every call too, and the complete event sequences, 164 to 1,541 lines
per case in the driver output, are equal for FeatureFinderCentroided_1 with its INI and with the
defaults (four charges), the four-scan input, a caller's map and a debug run
(`the_progress_event_sequence_matches_the_release_build`). They include the
inverted range the source passes in steps 2 and 3.2 on an input of fewer than
`2 * min_spectra` scans (`S 5 0`), which the port passes unchanged:
`ProgressLogger::start_progress` accepts it, as the Release build does
([PROGRESS_LOGGER_SUPPORT](PROGRESS_LOGGER_SUPPORT.md)).

## Reusing an instance

`instance::FeatureFinderAlgorithmPicked` is the source object with its state.
Its `run` is the source's `run`, and a second `run` of the same object behaves
as the source's second run (driver case reuse, three runs, and debug_twice):

- **The caller's map is extended.** `run` stores a pointer to the caller's map
  and clears it only for an empty input (`.cpp:134-138`, `1059-1063`). Step 3.3
  appends the new features, and step 4 sorts, resolves, filters and annotates
  the *whole* map: a caller's feature takes part in the overlap resolution
  (with the source's empty-box arithmetic for a hull without points and the
  full-plane box of a feature whose hulls include one), can become a
  subordinate or be removed, is re-sorted with the source's `std::sort`
  (`source_sort.rs`, a transcription of GCC 14's introsort, so ties land where
  the Release build puts them), and gets `spectrum_index` and
  `spectrum_native_id` from the current input. The map's unique id and meta
  values stay. Evidence: prefilled_run, prefilled_defaults_run and the five
  overlap cases, including NaN m/z, NaN intensity and infinite or NaN
  retention times, which the executed build completes.
  `algorithm::run` and `run_with_options` keep the fresh-map convenience.
- **`aborts_` accumulates**, and each run's abort block prints the
  accumulated counts.
- **`abort_reasons_` accumulates** and is keyed by intensity only
  (`Seed::operator<`), so seeds of equal intensity share one entry. A later
  debug run's abort map reads every stored seed's spectrum and peak index in
  *that* run's input (stale scaled: the first run's 25 seeds read on doubled
  intensities, 40 entries).
- **`isotope_distributions_` is extended, never cleared**
  (`IsotopeWindows::precalculate_onto`), so a reused object's windows, and with
  them its seeds and features, can differ from a fresh object's (reuse run 2).
- **`log_` is opened once and never closed.** A second debug run's `open`
  fails on the open stream, the stream's `failbit` drops every write, and the
  file keeps the first run's text (`DebugOutput::log_opened` is `false`). That
  holds after a first debug run that failed in steps 1 to 2.5 too: its `open`
  succeeded (`lenerr_reuse`, *Debug mode*). The first run's unflushed tail
  (up to 8,191 bytes) stays in the buffer until the object is destroyed; the
  instance keeps its counts (`debug_log_file`), and a process that ends in any
  later run, with or without `write_debug`, leaves the file at the first run's
  flushed length, which every termination reports (`log_file_bytes`; executed:
  *The file at termination*).
- **The user seeds are sorted in place.** `run_` sorts its member `seeds_`
  (`.cpp:190`), so after a run `seeds()` returns them in the Release build's
  introsort order; every run first replaces them with the caller's map.
- **Parameters.** `getParameters()` is the last run's merged set, or the
  last refused set: `setParameters`, which `run` calls first, assigns the new
  set before it checks it, and a refused set stays visible while the members
  keep the accepted values (`rejected_stdout.txt`: a direct call and a run),
  or, after a `ConversionError` in `updateMembers_`, with the members before
  the failing one updated (`partial_max_missing`, `partial_bins`).
  The warnings `checkDefaults` logs before it throws are the unknown keys it
  visits before the refused entry; a run puts them in its report.

Two undefined continuations of a reused object are refused at the point where
the source's process ends, and both record their termination: a stale abort
seed whose indices lie outside the current input (`Error::InvalidValue` while
building the abort map, `TerminationPoint::AbortMap`; the executed driver died
of SIGSEGV in all nine repetitions of rounds 1 and 5) and a caller's feature of
charge 0 in an overlapping pair with a different charge (`Error::InvalidValue`,
`TerminationKind::ArithmeticTrap`; the source's `%` traps, SIGFPE in every
repetition). The map sorts cannot read outside the map: their comparisons are
asymmetric, NaN keys included (*Non-finite input*, refusal 1).

## Preserved source conventions

- **Parameters.** The defaults carry the source names, values, descriptions,
  numeric and string restrictions and `advanced` tags. Their insertion order is
  the source's, and they equal the executed `-write_ini` section. Unknown
  entries are warnings, and type or restriction violations are errors, as in
  `DefaultParamHandler::setParameters`, with the source's texts and its
  32-bit narrowing of integers. `min_spectra` is
  `floor(mass_trace:min_spectra * 0.5)`, converted as the Release build
  converts it (the low 32 bits of `cvttsd2si`: `2^33 + 22` gives 11, `2^62 + 30`
  gives 0). The three percentages are divided by
  100. An abundance is changed when its `ParamValue` differs exactly from the
  default. Each changed abundance adds 1000 isotopes to the 20 of the patterns.
- **Check order of `run`.** The checks run in the source's order:
  1. no spectra returns an empty map before the parameters are read;
  2. no peak in spectra or chromatograms is an error;
  3. MS levels other than `{1}` are an error;
  4. unsorted spectra are sorted, with the source warning: the spectra by
     retention time and the chromatograms by product m/z with the Release
     build's `std::sort`, then each unsorted spectrum's and chromatogram's
     peaks with its `std::stable_sort` (`source_sort`). A spectrum (in
     retention-time order), then a chromatogram, that is sorted this way and
     holds a non-empty data array of another length than its peaks ends the
     run with the source's `Exception::Precondition` text, float arrays
     first, then string and integer arrays (`FloatDataArray[0] size (25) does
     not match spectrum size (24)`; executed: the `a_*` cases of
     `boundary_stage.tsv.gz`, ten of them with the text, four that run
     through);
  5. a negative first m/z is an error.

  The messages are the source's. The parameters are applied after these checks.
- **Intensity bins (step 1).** The grid spans the MS1 RT and m/z ranges.
  - The progress range is `startProgress(0, intensity_bins_ * intensity_bins_)`
    with a `UInt` product, which wraps modulo 2^32: 65,536 bins start `S 0 0`
    and 100,000 bins `S 0 1410065408` (executed, `progress_bins.tsv`); the
    `setProgress` values are `Size` and do not wrap (lead decision D12).
  - Bin borders are `start + i * step` and `start + (i + 1) * step`, evaluated
    in that form, and both are inclusive, as `areaBeginConst` is.
  - The area iterator visits only scans whose drift time lies in the full
    mobility range `[f64::MIN, f64::MAX]` (`MSExperiment.cpp:562-571`,
    `AreaIterator.h:277-298`): a NaN or infinite drift time leaves its scan
    out of every cell.
  - Intensities are promoted to `f64` and sorted with the Release build's
    `std::sort` (introsort), so signed zeros and NaN land where it puts them.
    Quantile `i` is element `floor(0.05 * i * (n - 1))`, and an empty cell
    keeps 21 zeros.
  - A peak's half-bin position is `floor((x - start) / step * 2)`, converted to
    `UInt` and capped at `2 * bins - 1`, and selects the neighbouring bins
    with the source's edge, odd and even rules. A zero or infinite step is
    computed, not refused (*Degenerate intensity bins*).
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
  - There are `ceil(max_mz * charge_high / width) + 1` windows, converted to
    `Size` as the Release build converts it; the progress starts with that
    count before the `resize`, which throws `std::length_error` above
    `vector::max_size()` = 164,703,072,086,692,425. Window `i`
    estimates the peptide mass `0.5 * width + i * width`, with at most 20 (or
    1020, 2020) isotopes, in binary32 source precision.
  - The source `trimLeft` keeps a pattern whole when no weight reaches the
    cutoff. It is followed by `trimRight`, both at
    `intensity_percentage_optional`; `trimmed_left` counts the removed leading
    isotopes.
  - Where all 20 binary32 bins of an estimate underflow, `renormalize`
    divides zero by zero and every weight is NaN: `trimLeft` erases nothing,
    `trimRight` erases everything, and the window is empty with maximum 0
    (executed from m/z `136,850.5` at charge 2 to `1e6`, and at charge 1000;
    the run continues). The executed boundary is an averagine mass between
    273,769.5 Da (the window centre at width 1 that still keeps one bin) and
    273,770.5 Da (the first empty one); at the FFC_1 width of 100 that makes
    window 2738 (centre 273,850 Da) the first empty one. The windows just
    below hold a single bin (`trimmed_left` 19): windows 2725 to 2737 at
    width 100 and 273,400.5 to 273,769.5 Da at width 1, as far as they were
    printed (`vw_w1`, `vw_w100` in `extended_stage.tsv.gz`). A NaN
    `intensity_percentage_optional`, which the parameter check accepts, fails
    the same two comparisons for every weight and empties every window
    (executed: `p_ipo_nan*`, no feature). The shared generator keeps its error
    for its other callers; step 2.5 uses the crate-private
    `estimate_from_peptide_weight_source`.
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
    - the correlation is `Math::pearsonCorrelationCoefficient` as the Release
      build computes it (`scoring::source_pearson`, crate-private): a
      denominator that underflows to zero from non-zero deviations gives an
      infinity, an all-equal range the default NaN;
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
    formed left to right, raised to `1.0f / 3.0f` by the reference build's
    `powf` (GNU C Library 2.39, `__powf_fma`), which is not correctly rounded
    everywhere: the port computes the same value, the misrounded ones and the
    NaN bits included (every binary32 base with this exponent executed). It is
    stored for every peak of the scored scans and compared with the thresholds
    in `f64`.
  - Seeds are local maxima at or above `seed:min_score`. With user seeds, the
    threshold is `user-seed:min_score`, and a peak also needs a user seed
    strictly within both user-seed tolerances. The search is the source's
    `lower_bound` plus forward walk over the seeds as the Release build's
    `std::sort` left them by m/z.
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
    `sqrt(max(0, 1 - mean relative deviation) * max(0, correlation))` and drops
    bad traces by the source's position rules relative to `max_trace`. The
    baseline is copied to the result last. The relative deviations are summed
    with the new term as the first operand, and every step follows the SSE
    operand order of `cropFeature_` and `checkFeatureQuality_`
    (`libOpenMS.so` `0x18d8a11`-`0x18d8ad4`, `0x18d9b5f`-`0x18d9e28`); the
    correlation is `source_pearson`'s.
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

## Degenerate intensity bins

Step 1 divides the MS1 retention-time and m/z extents by `intensity:bins`
(`FeatureFinderAlgorithmPicked.cpp:244-245`). The step is **zero** when every
MS1 spectrum has one retention time, when every MS1 peak has one m/z, or when a
subnormal extent underflows in the division (`4.9e-324 / 2`), and the
retention-time step is **infinite** when the extent overflows (retention times
from `-1e308` to `1e308`); either step is also infinite when a coordinate is
(see *Non-finite input*). The bins are still defined: with a zero step every
cell spans the extent, and with an infinite step the first cell's bounds are
`NaN` and `inf`, which the area iterator, like the port's inclusive walk, reads
as every spectrum. `intensityScore_` (`:1837-1838`) then converts
`floor(0 / 0)`, `floor(x / 0)` or `floor(inf / inf)` to `UInt` for every peak,
which is undefined behaviour.

**What the Release build does.** In `libOpenMS.so`
(`openms4-release-bc9cc12-c19e494-174b576`, GCC 14.4, `-O3 -mssse3
-ffp-contract=off`) the conversion is `cvttsd2si %xmm2,%rdi`, a 64-bit
truncation that gives `0x8000000000000000` for NaN and out-of-range values,
followed by the low 32 bits and an unsigned `cmovbe` cap. A NaN or infinite
position therefore selects half-bin 0; a position in `[-2^32, 0)` wraps to a
value above `2^31` and is capped at the last half-bin; a lower position keeps
its low 32 bits (at or below `-2^63`, half-bin 0). The cells read are always
inside the grid. The distances to the bin centres
are `0 / 0` or `inf / inf` whatever cell was selected, so every intensity score
is the default NaN `0xfff8000000000000` (stored as `0xffc00000`), every overall
score the seed loop computes is NaN, no seed is found, and the run returns an
empty map with the source's lines. Measured (three repetitions at one and at
four threads, identical):

| Evidence | Content | Test |
| --- | --- | --- |
| driver `degenerate_stage` (the B6 `seed_stage`, generalised), 26 configurations | FFC_1 with every RT equal, with every m/z equal, with a subnormal RT extent (1, 2 and 10 bins) and with an overflowing one, each with the FFC_1 INI, the defaults and `seed:min_score` 0; truncated to 10, 11, 14 and 15 scans, with and without an empty scan; the 499/500 m/z control; `FileFilter_44_input.mzML` and `FileConverter_31_output.mzML`. Members, bins, quantiles, windows, every float array of every peak with its NaN bits, `intensityScore_(spectrum, peak)` as a double, seeds, printed lines, feature counts and abort reasons | `degenerate_bin_steps_match_the_linux_release_build` |
| probe `iscore_probe`, four grids of 11 by 11 peaks | `intensityScore_` at 57 positions: in range, negative, beyond `2^31`, `2^32` and `2^63`, infinite, NaN, with NaN and infinite intensities, on a regular grid, a zero RT step, a zero m/z step and an underflowing RT step | `intensity_scores_outside_the_bins_match_the_linux_release_build` |
| FeatureFinderCentroided, cases `zero_rt`, `zero_mz` (each also at `-threads 4` and with `seed:min_score` 0), `zero_mz_control_min_score_0`, `filefilter_44_force` | exit 0, the printed lines, an empty map | `a_zero_width_retention_time_range_follows_the_cpp_release_build`, `a_zero_width_mz_range_follows_the_cpp_release_build`, `a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build` |

**What the port does.** `IntensityThresholds::score` follows the emitted
instructions: `x86_64::truncate_to_u32` is the conversion, and every arithmetic
step applies SSE2's NaN rule (the first NaN operand of the emitted instruction,
quieted, or the default NaN) in the operand order GCC emitted, which swaps some
commutative operands of the source text; the stored scores are narrowed with
the `cvtsd2ss` NaN rule. Non-NaN results are plain IEEE arithmetic and
unchanged. All 228 probe positions and every peak of the captures agree bit for
bit, NaN sign and payload included, on Linux x86_64 and on macOS arm64.

`Options::degenerate_bin_step` selects the behaviour:

- `DegenerateBinStep::Source`, the default and the tool's: compute as above.
- `DegenerateBinStep::Refuse`: `Error::InvalidValue` before any work, exactly
  when a step is zero or infinite **and** the seed loop visits a scan. The seed
  loop runs over `min_spectra_ .. n - min(min_spectra_, n)` (`:493-498`), which
  is empty for `n <= 2 * min_spectra_`; there the scores are never read and the
  result does not depend on them, so both variants return the source's empty
  map. Of the 26 captures, the 15 with a degenerate step and a non-empty seed
  loop are refused, and the other 11 run as under `Source`.

**The short input.** `FileFilter_44_input.mzML` has two MS1 spectra, both at
0.273 s, and the default `mass_trace:min_spectra` 10 makes `min_spectra_` 5, so
its seed loop is empty: its result is fixed by its length, not by its zero
extent. `FileConverter_31_output.mzML` (four spectra at 5 to 8 s, a non-zero
extent) behaves the same. The Debug oracle's exit 8 for both is
`OPENMS_PRECONDITION(begin <= end)` in `ProgressLogger::startProgress`
(`ProgressLogger.cpp:235`), reached from `startProgress(5, 0)` at `:298`, and
fires for every input shorter than `2 * min_spectra_` scans.

## Non-finite input

The source reads infinite and NaN retention times, m/z values, intensities,
drift times and user-seed positions without a check. The port follows it, and
refuses only where the source's behaviour has no reproducible answer (lead
decision D1: an out-of-bounds read or write, a data race, process termination
or a loop that never ends). What the executed Linux x86_64 Release build does
was measured with the driver `nonfinite_stage`
(`../oracle/ffap-sem-completion/drivers/nonfinite_stage.cpp`) on 189 modified
FeatureFinderCentroided_1 inputs, each run twice (identical; five also twice at
four threads, identical apart from the one-thread abort rows). The fixture
`nonfinite_stage.tsv.gz` holds every outcome, and
`non_finite_inputs_match_the_linux_release_build` replays all 189. The
combined fix round added 52 cases in the same format
(`sort_mobility_stage.tsv.gz`, `../oracle/ffap-complete-fix1/node/run_stage.sh`
with the drift-time variant `nonfinite_stage_dt.cpp`: every case twice,
identical, three also twice at four threads): the sorts over NaN, equal and
signed-zero keys, drift times and the step-2.5 bound, replayed by
`sort_and_mobility_cases_match_the_linux_release_build`. The third fix round
added 85 more (`boundary_stage.tsv.gz`, `../oracle/ffap-complete-fix3`, driver
`fix3_stage.cpp`, the same driver with options for swapped and reversed
spectra, data arrays and a chromatogram; every case twice, identical, two also
twice at four threads): underflowed averagine windows, a NaN trimming cutoff,
34 EGH and Gaussian fit configurations, the empty best pattern, the
charge-count wraps and mis-sized data arrays, replayed by
`boundary_cases_match_the_linux_release_build`. The fourth fix round
re-executed the round-3 numerics verifier's 124 cases and ten more
(`extended_stage.tsv.gz`, `../oracle/ffap-complete-fix4`, driver
`fix4_stage.cpp`, the verifier's `v3_stage`; every case twice, identical, and
the 124 shared rows equal to the verifier's own capture): fits on jittered,
skewed and rescaled inputs, isotope windows at 1-Da resolution across the
averagine underflow, both sides of every charge wrap and of the empty best
pattern, the `Precondition` order under ties, and the retention-time scales
whose `float` width overflows, replayed by
`extended_cases_match_the_linux_release_build`.

**What the source does, and what the port does.**

- *The input check* (`run`, `.cpp:1057-1098`). `MSExperiment::isSorted(true)`
  compares retention times with `>` and each spectrum with `std::is_sorted`;
  both are false for a NaN, so a NaN never makes the input unsorted. `-inf`
  sorts before every m/z and fails the positive-m/z check (`IllegalArgument`).
  The port runs the same comparisons (`validate_input`) and sorts as
  `sortSpectra` and `sortChromatograms` do, NaN and equal keys included
  (`source_sort`; executed: `v2_rt_*`, `v3_rt_*`, `v3_mz_nan_*`, and the
  2,272 inputs of `sort_probe.cpp`).
- *The ranges* (`MSExperiment::updateRanges`, `RangeBase::extend`). The
  minimum is `std::min` and the maximum `std::max` from `DBL_MAX` and
  `-DBL_MAX`: NaN never enters, the minimum is never `+inf` and the maximum
  never `-inf`. A range that stays empty (every retention time, or every m/z,
  NaN) throws `InvalidRange` ("Empty or uninitalized range object. Did you
  forget to call updateRanges()?") from `getMinRT`/`getMinMZ`. The port
  computes the ranges the same way and returns `Error::InvalidRange` with that
  text.
- *Step 1.* An infinite coordinate makes that step infinite, so every
  intensity score is NaN and no seed is found (*Degenerate intensity bins*;
  `DegenerateBinStep::Refuse` covers these inputs too). The area iterator
  finds its scans and peaks with libstdc++'s `lower_bound`/`upper_bound`; the
  port runs the same probe sequence (`scoring::libstdcxx`), so a NaN key gives
  the same cell contents. It skips a scan whose drift time is NaN or infinite
  (`v2_dt_*`, `v3_dt_*`: one scan, sixteen, all, and with an unsorted input;
  `f64::MAX`, `f64::MIN` and finite drift times keep their scan). Each cell is
  sorted with libstdc++'s introsort, NaN and signed zeros included (`v2_cell_*`,
  `v3_cell_*`). An infinite intensity sorts to the end of its cell and changes
  the quantiles; the score of that peak is `inf / inf`, NaN.
- *Step 2.* The trace scores compare intensities only, so an infinite
  intensity changes local maxima; a NaN or infinite m/z difference fails both
  comparisons of `positionScore_` (score 0). `findNearest` is reproduced with
  the same `lower_bound`.
- *Step 2.5.* `Size num_isotopes = ceil(max_mass / width) + 1` is converted by
  `comisd`/`cvttsd2si`/`btc` (`libOpenMS.so` `0x18e46f4`-`0x18e46fe`,
  `0x18e6c3b`-`0x18e6c44`): an infinite maximum m/z, or a count of `2^64` or
  more (m/z `1e300`), gives **no** window. A count above `vector::max_size()` =
  164,703,072,086,692,425 makes `resize` throw `std::length_error`
  (`vector::_M_default_append`; executed at `2^62 - 512`, `2^62`,
  164,703,072,086,692,448 and `1.5 * 2^63`), and the port returns that text.
  At or below it the source allocates `56 * count` bytes, which depends on
  memory (executed: 164,703,072,086,692,416 windows give `std::bad_alloc`);
  `Limits::max_isotope_windows` refuses every count above its ceiling there.
  The port converts with the crate-private `x86_64::truncate_to_u64`, which
  reproduces the sequence. A debug run has opened its log before either
  failure (*Debug mode*), and FeatureFinderCentroided reports the
  `length_error` from TOPPBase's `std::exception` handler, exit 12, as the
  executed tool does (`tool_1e19`); for the port's own ceiling it exits 8
  where the executed tool reports `std::bad_alloc` with exit 12 (`tool_2e18`,
  a native difference). Below the count limits, a window whose binary32 bins
  all underflow is emptied, as the source's NaN weights empty it (*Preserved
  source conventions*).
- *Step 3.1.* `getIsotopeDistribution_` (`0x18dddd0`) converts
  `floor(mass / width)` with the same sequence and throws `InvalidValue` when
  the index is not below the window count: with no window, at the first peak
  (`the value '12' ... Maximum allowed index is 0`; `'0'` for an infinite
  first peak), and for a NaN m/z at that peak, index `2^63`
  (`the value '9223372036854775808' ... Maximum allowed index is 15`). The port
  returns the same text (`IsotopeWindows::get`). With no charge
  (`charge_low = charge_high + 1`) step 3.1 never runs and the map is empty.
- *Step 3.2.* An infinite or NaN user-seed position never matches a peak
  (7 of the 8 FFC_1 features remain when seed 0 or seed 3 is changed); a
  single NaN seed matches nothing (no feature).
- *Step 3.3.* A NaN overall score is not below 0.01, so peaks with a NaN
  intensity score join mass traces; an infinite intensity makes the slope check
  cut the extension, and no executed feature holds a non-finite value. A NaN
  retention time in a mass trace makes `MassTraces::computeIntensityProfile`
  (`FeatureFinderAlgorithmPickedHelperStructs.cpp:210-236`) loop forever:
  every comparison is false, so the loop neither advances nor consumes a peak.
  The executed runs `rt_nan_mid` and `rt_nan_mid_bins3` did not return within
  30 s (killed, twice each; a stack sample shows the loop), which executes
  `CPP-242`. The port refuses at that merge (`MassTraces::intensity_profile`).
  With `feature:min_isotope_fit` 0 a seed whose best isotope pattern stayed
  empty is not aborted, and `extendMassTraces_` reads the first entry of the
  empty pattern (`pattern.spectrum[0]`): the executed build dies with SIGSEGV
  (`g_avg_trace0`, `g_iso0_seed0`, `p_ipo_100_seed0_iso0`,
  `p_ipo_nan_seed0_iso0`; a gdb backtrace shows the fault in
  `extendMassTraces_` under `run_`), and a bound of `1e-300` returns
  (`g_avg_trace0_iso_tiny`, 14 features). The port refuses at that read.

**Measured outcomes** (FFC_1 with its INI; `in`/`innear` sweeps set one peak of
or next to each of the 25 seeds):

| Cases | Executed Release outcome | Port |
| --- | --- | --- |
| `+inf`/`-inf` intensity at 3 fixed positions and 136 sweep positions (also with `seed:min_score` 0, `mass_trace:min_spectra` 2, the EGH model, `feature:min_isotope_fit` 0, reported m/z `average` and `maximum`) | 7 to 13 features; every score, seed, feature, meta value and hull recorded | the same (bit for bit on Linux x86_64; see below for macOS) |
| `-inf` m/z (first peak, last peak, a whole spectrum) | `IllegalArgument`, positive m/z | the same text |
| `+inf` m/z (last peak, a whole first spectrum, 5 spectra), m/z `1e300` | `InvalidValue`, no window | the same text |
| m/z `6.9e20`, `2.3e20` and below, `8.2e18` (window count above `vector::max_size()`) | `std::length_error`, `vector::_M_default_append` | the same text |
| m/z `8.2e18` with a window count just below `vector::max_size()` | `std::bad_alloc` | `Limits::max_isotope_windows` |
| `+inf` m/z with no charge | empty map | the same |
| `±inf` retention time (first, last, every; 3 bins; 10 scans) | empty map | the same |
| NaN m/z (one peak, a whole spectrum, a one-peak spectrum, a pair, both NaN signs) | `InvalidValue`, index `2^63` | the same text |
| NaN m/z in a whole spectrum with no charge | empty map, scores recorded | the same |
| every m/z or every retention time NaN | `InvalidRange` | the same text |
| NaN retention time of the first or last scan | 8 features | the same |
| NaN retention time of scan 50 (1 bin, 3 bins) | never returns | refused at the profile merge |
| user seed with `±inf` retention time, NaN retention time, `±inf` m/z | 7 features | the same |
| one user seed with NaN m/z; two, both NaN | no feature | the same |
| `mass_trace:min_spectra` 1, with and without an infinite intensity | empty map, trace scores `0xffc00000` (the default NaN of `0 / 0`) | the same bits on every host |
| NaN intensity (2 fixed positions, 6 sweep positions); NaN seed m/z among others | 8, 7, 7 features | the same, sorted as libstdc++ sorts them |
| NaN retention time with an unsorted scan (scan 50; scan 0) | never returns | refused at the profile merge |
| NaN or `-0.0`/`+0.0` intensities in one cell, strictly weakly ordered or not (12 to 40 one-peak scans) | no seed; quantiles and scores recorded | the same bits |
| equal retention times in an unsorted input (1, 3, 9 and 20 retention times set to another scan's); a NaN retention time at the second or last scan | 8, 8, 8, 7, 8, 8 features (`v2_rt_tie3_unsorted`: 26 seeds, only in the introsort order) | the same |
| unsorted spectra holding NaN m/z values | `InvalidValue`, index `2^63`; with no charge, the stable order in the trace scores | the same |
| user seeds with equal m/z, NaN m/z among equal and among different values | 0 to 5 features | the same |
| drift time NaN, `+inf` or `-inf` (one scan, 16, all), `f64::MAX`, `f64::MIN`, finite | 8 or 10 features, quantiles without the skipped scans | the same |
| last m/z `136,850` (2,738 windows) and `136,850.5` to `1e6` at charge 2 (the last 1 to 17,263 windows underflowed); also with `seed:min_score` 0, the EGH model, a zero optional cutoff, four threads; the first 20 scans at charge 1000 | 8, 13 and 0 features | the same, the underflowed windows empty |
| `intensity_percentage_optional` NaN (FFC_1, `seed:min_score` 0, user seeds with both thresholds, 3 bins, the EGH model), 100 and `-0.0` | 0 features; 8 at `-0.0` | the same, every window empty |
| `feature:min_isotope_fit` 0 with `seed:min_score` 0 (and reported m/z `average` with the other thresholds 0; an optional cutoff of 100 or NaN) | SIGSEGV, twice each | refused at the empty best pattern |
| `charge_low`/`charge_high` 3/2, 4/2, 5/2, 6/2, 7/2, `INT_MAX`/1, `INT_MAX`/498, 1/`INT_MAX`, 2/`INT_MAX` | empty map; SIGSEGV; `std::bad_alloc` three times; SIGSEGV three times; `std::bad_alloc` | empty map; refused as undefined; refused before the allocation; refused as undefined; the native charge ceiling |
| unsorted input with a float, string or integer data array of another length than its peaks, in one or two spectra (also with every retention time equal), or in a chromatogram; sorted spectra and exact arrays | `Exception::Precondition` for the first such spectrum in introsort order, then chromatograms; otherwise 8 features (0 with one retention time) | the same text; the same |

**What is refused, and where.** Each refusal is `Error::InvalidValue`, at the
first point where the source's behaviour has no reproducible answer:

1. *A step of `std::sort` that would read outside the vector: an unreachable
   guard here.* libstdc++'s introsort (`bits/stl_algo.h`) has no bound checks
   in its partition and final insertion loops. The port runs the same
   comparisons and moves and would refuse exactly at a read outside the
   vector (`source_sort`), but every sort of this algorithm compares with an
   asymmetric `<` (on `f32` or `f64` keys, NaN included), and for any
   deterministic asymmetric comparison both loops provably stay inside the
   vector (the argument is in the `source_sort` module documentation). No
   executed or generated input reaches the guard: not the 189 and 52 stage
   cases, not the 130 inputs of `ffap_instr_driver` (mode `sort`), not the
   2,272 of `sort_probe.cpp`, not the 10,275 cases of the numerics verifier's
   probe, and not the unit test's 4,000 NaN key sets and 1,000 asymmetric
   tournaments. Only a comparator that is not asymmetric, such as `<=`, reaches
   it through the public `source_sort_permutation`; the port then follows the
   source's signed iterator arithmetic while it reads inside the vector and
   never panics (`other_comparisons_are_refused_or_sorted_without_a_panic`).
   Everything else about a NaN or equal sort key is reproduced: the step-1
   cells (`.cpp:270`), the spectra and chromatograms (`sortSpectra`,
   `sortChromatograms`), the user seeds (`.cpp:190`), the seeds (`.cpp:548`)
   and the feature map (`.cpp:866`, `:991`). `std::stable_sort` never reads
   outside its range, whatever the comparisons return, so it is reproduced on
   every key: its merges check both ends, its binary searches and rotations
   are bounded, and its insertion sort stops at the first element, which the
   inserted one was just found not to be less than.
2. *An area-iterator cell whose lower search lies above its upper search*:
   `AreaIterator` would never meet its end and read past the scans or the
   peaks (`AreaIterator.h:205-218`, `:276-298`). This cannot occur, even with
   NaN keys: a cell's lower border never exceeds its upper one (or one is
   NaN), and for `v <= w` libstdc++'s `lower_bound(v)` and `upper_bound(w)`
   branch alike at every probe until the first probe where they differ, which
   sends the lower search left of it and the upper search right of it. The
   port keeps the check as the guard of its slices; no input reaches it.
3. *A NaN retention time merged into an intensity profile*: the source never
   terminates (above; `CPP-242`).
4. *An empty best isotope pattern read at `extendMassTraces_`*
   (`.cpp:1347-1349`): with `feature:min_isotope_fit` 0 a seed without any
   placement reaches the read of `pattern.spectrum[0]` of an empty vector, an
   out-of-bounds read; the executed build dies with SIGSEGV (four stage cases,
   the three negative-intensity library runs of *Debug mode* with and without
   `write_debug`, and FeatureFinderCentroided). Reachable with ordinary
   parameter values. A debug run keeps the seed's lines and records a
   `TerminationKind::OutOfBounds` termination.
5. *The `UInt` score-array count that wraps* (`.cpp:196-221`, lead decision
   D12): refused for every wrapping count whatever the `Limits`, as
   `Settings::charge_count` documents; executed at seven wrapping pairs. The
   counts `-2` and `-3` stay refused unconditionally (lead decision D13):
   their executed outcome, `std::bad_alloc`, depends on memory. Where the
   process dies there (an out-of-bounds write after an allocation the port
   takes to succeed), the run records a `ScoreArrays` termination, before
   `debug/log.txt` is opened (*Debug mode*).
6. *A feature m/z without an isotope window at step 3.3.5* (`.cpp:790`): the
   source's `InvalidValue` leaves its OpenMP region uncaught and
   `std::terminate` ends the process. A NaN feature m/z (an infinite intensity
   kept in the reported traces), or an `average` m/z whose intensity sum
   nearly cancels, would reach it; no executed input did (source review; the
   instrumentation verifier's 50 generated inputs with infinite isotope
   intensities all returned, and fix round 5 searched the port, which
   reproduces the executed runs, over 18,720 FeatureFinderCentroided_1
   variants with one m/z band of intensities negated by factors from 0.05 to
   40, with `reported_mz` `average` and `maximum`, `feature:min_isotope_fit`
   `1e-300`, the other thresholds 0 and `seed:min_score` 0 and 0.3, and 35,100
   coarser variants before: none reached `.cpp:790`, so no candidate was
   executed). In a debug run the source has written that seed's log lines and
   feature files first, and the port keeps them before it records the
   termination.

**No platform split.** Since lead decision D10 both trace fitters call the
reference build's glibc `exp` and `log`, ported (`glibc_libm`), so every
returned feature of the three stage fixtures is bit for bit on Linux x86_64
and on macOS arm64, the ill-conditioned fits that used to depart on macOS by
up to `4.9245e-4` included. The one host-dependent call left is the `atan` of
an EGH area (`glibc_libm::atan`, D10's fallback): on a host other than x86_64
Linux with the GNU C Library the port calls the `libm` crate's, and over these
fixtures that moved no EGH intensity on macOS arm64 (a measured maximum of 0,
not a guarantee; `area_tolerance` in the test; rates under *Platform
notes*).

**The tool.** `FeatureFinderCentroided` reads the mzML with the native reader,
which refuses non-finite binary values and scan start times; see
[TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md),
native difference 13.

## Native differences

| Source behaviour | This port | Reason and evidence |
| --- | --- | --- |
| scores are float data arrays appended to each spectrum, replacing its existing float arrays | `ScoreArrays`: one flat `f32` array per score, outside the spectra | Only the algorithm reads them, and only debug mode writes them. The input spectra keep their arrays, and no per-spectrum allocation is needed. |
| `spectrumRanges().byMSLevel(1)` needs `updateRanges()` from the caller | ranges are computed on demand from the validated spectra, with `RangeBase`'s `std::min`/`std::max` semantics | No stale-range state exists, so the source's FAIMS "No ranges for this MS level" crash cannot occur. The source message "needs updated ranges" belongs to the peak-count check and is kept verbatim. |
| infinite and NaN retention times, m/z values, intensities, drift times and user-seed positions are read without a check | read as the Linux x86_64 Release build reads them; refused only at the endless profile merge of a NaN retention time and at the uncaught step-3.3.5 exception (the introsort's out-of-bounds guard is unreachable for these keys) | See *Non-finite input*: of the 189 executed cases 170 returned runs and 16 exceptions are reproduced and the 3 endless runs refused at the merge; of the 52 fix-round cases 43 returned runs and 8 exceptions are reproduced and 1 endless run refused. |
| `mass_trace:min_spectra = 1` gives `min_spectra_ = 0`. Every trace score becomes 0/0 = NaN and every peak a local maximum; no overall score reaches a threshold, no seed is found and the run returns an empty map | the same: NaN trace scores, no seed, an empty map | The execution (B6 driver, `ffc1_min_spectra_1`) shows the source is *defined* here, so the port follows it (lead decision of 2026-09-15, `CPP-271`). B6 refused the configuration; that refusal is gone. Nothing later in the algorithm is reached, so the source's `size_t(-1)` delta buffer in `extendMassTrace_` stays unreachable; the port returns `Error::InvalidValue` if it ever is. |
| a zero or infinite intensity bin step makes `intensityScore_` convert `floor(NaN)` or `floor(inf)` to `UInt`, which is undefined | by default the Linux x86_64 Release build's outcome (every intensity score NaN, no seed, an empty map); `DegenerateBinStep::Refuse` refuses exactly the inputs whose seed loop reads those scores | Undefined behaviour of the `float`-to-`int` kind, whose Release outcome is measured, repeatable and explained by the emitted `cvttsd2si`; see *Degenerate intensity bins*. The port refused every zero-width range before, including short inputs, where the result does not depend on the scores. |
| `UInt charge_count = charge_high - charge_low + 1` and the `UInt` array count `3 + 2 * charge_count` wrap for `charge_low > charge_high + 1` and for `charge_low` 1 with `charge_high` `INT_MAX` | `Error::InvalidValue` from `Settings::charge_count` for every wrapping count, whatever the `Limits` | Lead decision D12. The counts `2^31 - 1`, `-1` and `-4` and below wrap to an array count the pattern loop writes past: an out-of-bounds write, undefined (executed SIGSEGV: 4/2, `INT_MAX`/1, `INT_MAX`/498, 1/`INT_MAX`). The counts `-2` and `-3` wrap to `2^32 - 1` and `2^32 - 3` arrays per spectrum, which stay in bounds if the allocation succeeds; that outcome depends on memory (executed `std::bad_alloc`: 5/2, 6/2), so the port refuses before allocating, and before a debug run opens its log, as the executed runs had not opened it either. `charge_low == charge_high + 1` gives zero charges, as in the source (executed: an empty map). |
| `write_debug = true` writes `debug/` into the working directory, and `writeFeatureDebugInfo_` reads the undeclared `debug:pseudo_rt_shift`, which terminates the process at the first seed that reaches the fit | `debug::DebugOutput`, which the tool writes; `Error::Unsupported` at that seed under `PseudoRtShiftKey::Source` | Library code does not write files, and a safe port cannot end the process (see *Debug mode*). |
| a changed `abundance_12C` or `abundance_14N` builds the override from a default `IsotopeDistribution` that already holds `(0, 1)`; the patterns grow (FFC_1 window 0: 27 normalised bins) and FFC_1 with 12C = 90 % finds 0 seeds, 0 candidates and 0 features (C2 `ffap_ffc1_abundance_12C_90`) | the intended two-isotope distribution, which **does** find seeds and features | `CPP-247`, lead decision of 2026-09-15: follow the intent, not the defect, because the generator rejects the stray-peak construction and refusing a parameter the source accepts is worse. **This is the one place where the port's features differ from the executed C++ by design.** `AbundanceOverride::Refuse` is the opt-in for a caller that must not diverge. |
| `std::sort` of seeds with equal `f32` intensity leaves an order the standard does not specify | the order the Linux Release build's libstdc++ introsort gives (`source_sort::source_sort_reversed_by`) | Deterministic, and equal to the reference platform (driver mode `sort`, 130 inputs). |
| the overall score is `std::pow(float, float)`, the platform's `powf` | the reference build's `powf`, ported (`glibc_powf`, crate-private), on every host | Lead decision D5: the Linux x86_64 Release build binds `powf@GLIBC_2.27` of glibc 2.39, which selects `__powf_fma` on the reference CPU; the port reproduces its instructions, the fused multiply-adds included, and the results equal the executed library for every binary32 base with the exponent `1.0f / 3.0f` (2^32 values, 256 digests), for `2^26` generated pairs in each of four sets and for a 42 by 42 special grid (`../oracle/ffap-complete-fix1`, `powf_probe`; the unit tests replay the grid, 1,024 rows and `2^20`-pair digests). All 30,840 retained overall scores and every score of the stage fixtures are the executed ones; the 8 retained scores (2 of 3,084 on FFC_1) that glibc rounds one binary32 step from the correctly rounded power are listed in `overall_rounding.tsv` and asserted as misrounded. The macOS arm64 product SDK's Apple `powf` misrounded 99 others; the port no longer depends on the host's `powf`. |
| progress, `Found N seeds for charge c.` and `Found N feature candidates for charge c.` go to `std::cout`; the overlap count, the abort reasons, the feature count and the apex warning to `OPENMS_LOG_INFO`/`WARN` | `SeedStage::log` and `RunOutput::log`, in the source's order | Library code never prints. The candidate line is inserted directly after its charge's seed line, so the two `std::cout` lines are adjacent as in the source, even though this port computes every charge's seeds first. The bare newlines the source logs around the abort block (`FeatureFinderAlgorithmPicked.cpp:1019` and `1026`) are emitted as empty log entries, so a caller that prints the log line by line — `FeatureFinderCentroided` does — reproduces the executed C++ stdout block exactly. |
| steps 3.1 to 3.3 run per charge | steps 3.1 and 3.2 run for every charge first | Step 3.3 reads only its own charge's arrays, so the arrays and seeds are identical. |
| user seeds are a copied `FeatureMap` sorted in place with `std::sort` | the same map sorted in place in the Release build's introsort order; the seed stage keeps positions only | The m/z and retention time are the only fields the source reads. Equal and NaN m/z values are sorted as the executed library sorts them (`sort_probe.cpp` site `features`, and the `v2_seeds_*`/`v3_seeds_*` runs). |
| unbounded work | `Limits`: spectra, peaks, charges, bins per dimension, windows, pattern values, score bytes, work units | Checked before the allocation or computation each bounds. The FFC_1 workload is several orders of magnitude below every default. |
| `Math::pearsonCorrelationCoefficient` returns an infinity when its denominator underflows to 0 from non-zero deviations | the same infinity (`scoring::source_pearson`, crate-private, for the isotope, crop and quality correlations) | Since fix round 3; the crate's shared `pearson_correlation_coefficient` keeps its NaN for its other callers. In `cropFeature_` and `checkFeatureQuality_` `std::max(0.0, inf)` keeps the infinity, and the final score is infinite or NaN (`inf * 0`), which no threshold rejects. Reaching it needs theoretical intensities that vary by less than about `1e-154` around the baseline while the measured ones vary; no executed or generated input did. Such a feature then stores the non-finite `score_correlation`, as the source does (next rows). |
| `Exception::UnableToFit` from `fitter->fit` would be thrown inside the `omp parallel for`, where nothing catches it, and would end the process | unreachable, so not reproduced; an error from the port's fit (its point, byte or work ceilings, or the refused NaN profile merge) ends the run with that error | Both throws of `TraceFitter::optimize_` are unreachable from the seed loop: a fitted candidate has at least two traces of which at most one has fewer than three peaks, so at least 4 residuals for at most 4 parameters (`TraceFitter.cpp:111`), and Eigen returns `ImproperInputParameters` (`:129`) only for `maxfev <= 0`, which the `fit:max_iterations` restriction excludes. The residual count is `int` (`GaussTraceFitter.cpp:140`, `EGHTraceFitter.cpp:29`), so more than `INT_MAX` peaks would also throw at `:111`, but the solver's `MAX_POINTS` ceiling refuses such traces first. The argument is written out at `FittedModel::fit`. The port used to turn any fit error into an abort reason, which hid its own ceilings. |
| `extendMassTraces_` reads `pattern.spectrum[0]` when the best pattern stayed empty (or its first matched isotope has no peak) | `Error::InvalidValue` at that seed (`extension::EMPTY_PATTERN_WHAT`, crate-private), with a `TerminationKind::OutOfBounds` termination in a debug run | Undefined behaviour: an out-of-bounds read (lead decision D1). The SIGSEGV and the flushed-only log are established for an **empty** best pattern, which every executed crash had. The other sub-case, a non-empty pattern whose first isotope has no peak, reads `map_[spectrum][size_t(-1)]`, 16 or 32 bytes before a spectrum's peak buffer (heap metadata, which is mapped), so the process would probably not fault there; it was never observed (the round-3 instrumentation verifier classified 6,076 refusals of a 19,200-run grid: all empty patterns), and its executed outcome is unknown. The port refuses it at the same read and records the same termination, whose `SIGSEGV` label is established only for the empty pattern. Reachable from `run_` with ordinary parameter values: `feature:min_isotope_fit` 0 lets a seed without any placement through. Executed SIGSEGV, twice each: the stage cases `g_avg_trace0`, `g_iso0_seed0`, `p_ipo_100_seed0_iso0` and `p_ipo_nan_seed0_iso0`; the library runs `neg_oob1`, `neg_oob_seed035` and `neg_none_avg0` (negated intensities or none, with and without `write_debug`); FeatureFinderCentroided `avg0` (`../oracle/ffap-complete-fix3`; `../oracle/ffc-numerics-v2`, `logs/g_avg_trace0_gdb.txt`, and `../oracle/ffc-instrumentation-v2` found the first ones). `boundary_cases_match_the_linux_release_build`, `a_seed_loop_crash_keeps_what_the_executed_process_had_written`, `a_run_that_reaches_an_empty_best_pattern_is_refused_where_the_release_build_crashes`; the public function's refusal: `an_empty_pattern_is_refused_instead_of_dereferenced`. |
| `traces[traces.max_trace]` is indexed before the fit without a range check | `Error::InvalidValue` | Undefined behaviour when `max_trace` is stale. The one branch that could make it stale (`traces.clear()` for a trace before `max_trace`) is unreachable, because `max_trace` is still 0 at every index that could satisfy `p < max_trace`. |
| `setWidth` stores any FWHM, and `setMetaValue` any `score_fit`, `score_correlation` or `EGH_*` value | the same values, non-finite and negative ones included: the `width` field directly, the meta values through the crate-private `MetaValue::source_float` (accepted by lead decision D13); `BaseFeature::validate`, `MetaValue::validate` and the featureXML writer refuse them | Since fix round 4 (lead decision D1: measured, repeatable, explained and in bounds). Reachable with finite input and unchanged FFC_1 parameters: once the fitted `sigma` passes about `1.44e38` the `float` FWHM overflows, and the Release build returns features with an infinite width, `FWHM` meta value and intensity. Executed onset on FeatureFinderCentroided_1 with every retention time scaled (of 9 Gaussian and 8 EGH features): no infinite width at `2e36` and `4e36`, 1 and 0 at `6e36`, 3 and 2 at `8e36`, 7 and 6 at `1e37`, all from `1.5e37` (`width_onset_stage.tsv.gz`, `extended_stage.tsv.gz`: `vy_rt_1e37` to `vy_rt_1e150`, `vx_rt_1e150` to `vx_rt_1e300`); every intensity is infinite from `1e36` on and finite at `1e33`. The round-3 numerics verifier found it; before, the port refused the run. The public `TryFrom<f64>` of `MetaValue` still rejects non-finite values; only the algorithm's own stores skip the check. A non-finite `score_correlation` (the zero denominator above) and a NaN FWHM are stored the same way; neither was reached. The debug seed map stores its scores through the same constructor only to follow `setMetaValue(float)`'s storage: a seed's scores are always finite (`debug::seed_map`). FeatureFinderCentroided cannot write a map with a non-finite feature value (TOPP native difference 16; the featureXML writer's refusal and the CLI's wording for that write failure are split off into a separate task, lead decision D13). |
| `f2.getCharge() % f1.getCharge()` divides by zero for a zero charge | `Error::InvalidValue` at that pair, only where the remainder the source evaluates would trap (also `INT_MIN % -1`) | The algorithm's own charges are at least 1. A caller's feature of charge 0 traps in the source (executed: SIGFPE, 2 of 2); same-charge pairs and non-overlapping charge-0 features are processed (executed). |
| a hull with no point has the default `DBoundingBox` `[DBL_MAX, -DBL_MAX]`, whose `width()` is negative infinity | the same arithmetic (`resolution.rs` `SourceBox`) | Only a caller's map can hold such a hull (executed: overlap cases). |
| `plot_nr` is assigned in an OpenMP critical section, so its value depends on the schedule | assigned in seed order | It is overwritten by the feature number for every feature that survives. The debug file names use it, and with one thread the source's value is this seed-order number. |
| `aborts_[reason]++`, `abort_reasons_[seed]` and the `log_` writes of `abort_` run inside the parallel region without synchronisation | aggregated serially in seed order | A data race (candidate 5; `.cpp:595`, calls at `:627`, `:640`, `:725`; the executed tool at four threads wrote a different log in each run, cases c1 and c2): with two or more threads aborting seeds the source has no reproducible result. **Lead decision D11 accepts this as the one documented exception to D1**, which lists data races as a refusal class: the port gives the source's single-thread result at every thread count, because the user's determinism contract requires parallel output to equal serial output, and refusing would block essentially every parallel run, as `aborts_` is written in every run. Every executed comparison of `aborts_` (C2, `degenerate_stage`) records the library's map at one thread only, and the multi-thread runs compare everything else. |
| the seed loop is an OpenMP `parallel for` with four named critical sections | `concept::parallel::map_collect` over the seed indices with `Options::threads` | `map_collect` returns results in input order, so no critical section is needed and the output is bit-identical at 1, 2 and 8 threads (`the_seed_loop_is_bit_identical_across_thread_counts`). The source's results are schedule-independent for the same reason, its `tmp_feature_map` being a `std::map` keyed by seed index; the C++ oracle gave identical output at 1, 2 and 8 threads. |
| the containment pass scans a growing `std::vector<Size>` of swallowed seeds | a `BTreeSet` | A pure membership test; duplicates in the source's vector change nothing. |
| `FeatureMap::sortByMZ` and `sortByIntensity` use `std::sort`, which is unstable | the Linux Release build's introsort order (`source_sort::source_sort_by`), NaN keys included | Equal to the reference platform (driver modes `sort`, `overlap`); the out-of-bounds guard is unreachable for these comparisons. |
| `setMetaValue("spectrum_index", Size)` stores an unsigned value | `i64`, refused above `i64::MAX` | The source's `DataValue` narrows to a signed integer anyway; the featureXML writer writes the same digits. |
| `CoarseIsotopePatternGenerator` iterates elements in heap-address order, which varies between runs | ascending atomic number | Inherited from B2. It is the majority order (198 of 200 runs), which every retained execution used; all windows match. |
| `MSExperiment::sortSpectra` and `sortChromatograms` | module-local sorts that reproduce them: libstdc++'s introsort for the spectra and chromatograms, and per unsorted spectrum or chromatogram its `std::stable_sort` (`source_sort`), NaN and equal keys included | The kernel's sorts refuse every non-finite value. A spectrum or chromatogram whose data arrays do not match its peaks gives the source's `Exception::Precondition` text at the source's point, the first such unsorted spectrum in introsort order (executed: ten `a_*` cases). The spectra sorted before it stay sorted in the source's caller-visible map; `run` consumes the port's experiment, so no caller can look at it. FeatureFinderCentroided cannot reach it: the native mzML reader refuses an array length other than `defaultArrayLength`. |
| `std::stable_sort` asks `operator new(nothrow)` for its buffer and halves the request after each failure | the same halving on the port's own allocation (`TemporaryBuffer::Allocate`); the port's buffer holds 8-byte indices where the source's holds 16-byte peaks or 8-byte indices | Which requests fail depends on the process's memory and is not reproducible; the algorithm for every buffer size is the source's (executed with a replaced `operator new(nothrow)` that refuses above a byte limit: full, partial and no buffer, `sort_probe.cpp`). |
| step 2.5 allocates `56 * count` bytes for a window count at or below `vector::max_size()` | `Limits::max_isotope_windows` (default 1,000,000) refuses larger counts first, after the point where a debug run has opened its log; FeatureFinderCentroided exits 8 with that message where the executed tool reports `std::bad_alloc` with exit 12 | Allocation failure depends on memory (executed: 164,703,072,086,692,416 windows, and the `4e16 + 1` of m/z `2e18` at charge 2 and width 100, throw `std::bad_alloc` on the reference node); above `max_size()` the port returns the source's `length_error` text and the tool exits 12 as the executed one does (lead decision D6; m/z `1e19` gives `2e17 + 1` windows). |
| the charge loop resizes every spectrum's float arrays to `3 + 2 * charge_count` in `UInt`, `charge_count` up to `2^31 - 2` without wrapping | `Limits::max_charges` (default 1,000) refuses more charges | A native ceiling in front of an allocation that depends on memory (executed: `std::bad_alloc` for `charge_low` 2, `charge_high` `INT_MAX`). The wrapping counts are refused whatever the limits (*Non-finite input*, refusal 5). |
| `MSExperiment::RTBegin`, `MSSpectrum::MZBegin`/`MZEnd`/`findNearest` and the quantile search use `std::lower_bound`/`upper_bound` | the same probe sequence (`scoring::libstdcxx`) | On sorted keys any binary search gives the same index; on keys a NaN leaves unpartitioned only libstdc++'s sequence does. |

### Score arrays

The source names are kept by `ScoreArrays::array_names`: `trace_score`,
`intensity_score`, `local_max`, `pattern_score_<c>` for each charge, then
`overall_score_<c>` for each charge. Trace scores, local-maximum flags and overall
scores stay zero outside `min_spectra..len - min_spectra`, as in the source.

## Checked boundaries and evidence

### Tier 1: executed C++ (Linux x86_64 Release)

The reference platform is the Linux x86_64 Release build
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576` (core
`bc9cc12`, cli `c19e494`, topp `174b576`; `libOpenMS.so` sha256 `abd4fc99...`),
the user's decision of 2026-09-15. The seed-stage and feature-stage fixtures
were first captured from the macOS arm64 product SDK (Debug, core `4fdec46`,
fused-FMA Eigen, Apple libm) and were re-captured from the Release build with
the same drivers on ibminode06 (AMD EPYC 7763, glibc 2.39); the B6 and B7
extraction scripts ran unchanged on that capture
(`../oracle/ffap-sem-completion/extract`). Only the `-write_ini` comparison
still uses the product SDK's executed file.

Every float of the **seed stage**, the overall scores included, is compared bit
for bit on every platform. The **feature stage** compares counts,
identities and orders exactly, and fitted quantities bit for bit on every
platform, except the area and intensity of an EGH fit on a host other than
x86_64 Linux with the GNU C Library (`area_tolerance`, a measured maximum of
0); the measured agreement is in the next section.

| Evidence | Configurations | Test |
| --- | --- | --- |
| C1 `FFC_write_ini`: FeatureFinderCentroided `-write_ini` (product SDK, two runs, identical) | the whole algorithm section | `default_parameters_equal_the_executed_write_ini_section` |
| C2 `ffap_stages` library state (`omp1` and `omp4`, two repetitions each, identical): effective members, bins, 21 quantiles, 15 windows, all score arrays of 3,084 peaks, printed seed counts | FFC_1 INI (tool loading), class-test INI, #9247 tight pattern and tight trace, FFC_1 with the retained output as user seeds | `ffc_1_stage_...`, `class_test_stage_...`, `tolerance_swap_stages_...`, `user_seed_stage_matches_the_executed_library` |
| B6 `seed_stage` library state, as the generalised `degenerate_stage` (three repetitions, identical after dropping its extra row kinds), same content | default parameters (10 bins, charges 1 to 4, width 25, 107 windows); FFC_1 INI with 7 bins and charges 1 to 3 (21 windows) | `default_parameter_stage_...`, `seven_bin_three_charge_stage_matches_the_executed_library` |
| the same, `ffc1_min_spectra_1` | `mass_trace:min_spectra = 1`: `min_spectra_` 0, 0 seeds, 0 features, exit 0 | `min_spectra_one_follows_the_source_and_finds_no_seed` |
| C2 `ffap_stages` **final `FeatureMap`**, `aborts_` and the two `std::cout` lines, for six configurations | FFC_1 symmetric, FFC_1 asymmetric (EGH), FFC_1 with the retained output as user seeds, the class-test INI, the two `#9247` tolerance swaps | `every_configuration_matches_the_executed_library` |
| `degenerate_stage` and `iscore_probe` | the degenerate intensity bins (26 configurations, 228 probe positions) | see *Degenerate intensity bins* |
| `defs_probe` (three repetitions, identical) | `FeatureFinderDefs` | `feature_finder_defs_match_the_executed_probe` |
| `nonfinite_stage` (two repetitions each, identical; five cases also at four threads) | 189 FFC_1 inputs with infinite or NaN retention times, m/z values, intensities and user-seed positions | `non_finite_inputs_match_the_linux_release_build`; see *Non-finite input* |
| `nonfinite_stage_dt` (`../oracle/ffap-complete-fix1`, two repetitions each, identical; three cases also at four threads) | 52 FFC_1 inputs: sort keys (NaN, equal, signed zero), drift times, window counts around `vector::max_size()` | `sort_and_mobility_cases_match_the_linux_release_build` |
| `sort_probe` (two runs, identical) | `MSSpectrum::sortByPosition` and `MSChromatogram::sortByPosition` with and without a data array under a full, a partial and no temporary buffer; `sortSpectra`, `sortChromatograms`, `FeatureMap::sortByMZ`; 2,272 inputs | `every_source_sort_matches_the_executed_library` |
| `powf_probe` (two runs, identical) | the resolved `powf`, a 42 by 42 special grid, four generated sets, every binary32 base with the exponent `1/3` | unit tests of `glibc_powf`; the exhaustive and `2^26`-pair digests were compared outside CI with a release build of the same module (`../oracle/ffap-complete-fix1/port-harness`) |
| `defs_eq_probe` (two runs, identical) | `ChargedIndexSet` comparisons | `charged_index_set_comparisons_match_the_executed_probe` |
| `fix2_driver bigint` (`../oracle/ffap-complete-fix2`, two runs, identical) | 64-bit values of `intensity:bins`, `isotopic_pattern:charge_low`/`charge_high`, `mass_trace:max_missing`/`min_spectra` and `fit:max_iterations` beyond the `int` range: outcomes with their texts, features bit for bit, aborts, console lines, the protected members, the stored values; two refused sets after an accepted one; a negative `fit:max_iterations` at the start of `run_` | `integer_parameters_beyond_the_int_range_follow_the_release_build` |
| `fix2_driver lenerr_single`, `lenerr_reuse` and the Release FeatureFinderCentroided (two runs each, identical apart from timing text) | m/z `1e19` (`length_error`) and `2e18` (`bad_alloc`) in a debug run: the exception, the debug directory and log while alive and after destruction, a reused object's second run against a fresh object's; the tool's exit, stderr, stdout block and files | `a_debug_run_that_fails_in_step_two_point_five_keeps_the_executed_debug_output`, `a_debug_run_beyond_the_isotope_window_limit_exits_as_the_release_build` |
| `fix2_driver formula` (two runs, identical) | `GaussTraceFitter` and `EGHTraceFitter::getGnuplotFormula` over 12 special values (signed zeros, infinities, four NaN bit patterns) in the sum and the product, and in `sigma`, `tau` and the baseline: 648 formulas | `gnuplot_formulas_print_the_executed_nan_signs` |
| `libm_probe` (`../oracle/ffap-complete-fix3`, two runs, identical) | the `exp`, `log` and `atan` `libOpenMS.so` binds (glibc 2.39: `__ieee754_exp_fma`, `__ieee754_log_fma`, `__atan_fma`): an 80-value special grid, and 13 generated sets of `2^26` inputs (64 digests and 256 rows each) | unit tests of `glibc_libm` (the grid, the rows and the first `2^20` of each set; `atan` only with glibc); every input outside CI with a release build of the same module (`port-harness`: `exp` and `log` equal on macOS arm64 and on kim, `atan` equal on kim) |
| `fix3_stage` (two runs each, identical; two cases also at four threads) | 85 FFC_1 inputs: underflowed averagine windows, a NaN trimming cutoff, 34 EGH and Gaussian configurations, the empty best pattern (SIGSEGV), the charge-count wraps (SIGSEGV, `std::bad_alloc`), unsorted input with mis-sized data arrays (`Exception::Precondition`) | `boundary_cases_match_the_linux_release_build`, `charge_count_wraps_are_refused_whatever_the_limits` |
| `vfi2_driver neg` (after `../oracle/ffc-instrumentation-v2`; two runs each, identical) | three negative-intensity inputs with `feature:min_isotope_fit` 0, with and without `write_debug`: SIGSEGV; the debug log and every feature file left behind | `a_seed_loop_crash_keeps_what_the_executed_process_had_written` |
| `fix3_driver progress` (two runs each, identical) | the step-1 `startProgress` event for nine `intensity:bins` values, 65,536 and 2^32 + 65,536 among them | `the_step_one_progress_range_wraps_as_the_release_build_computes_it` |
| the Release FeatureFinderCentroided, case `avg0` (two runs, identical apart from timing text) | the empty best pattern through the tool: SIGSEGV, the console lines written before it | `a_run_that_reaches_an_empty_best_pattern_is_refused_where_the_release_build_crashes` |
| `fix4_stage` (`../oracle/ffap-complete-fix4`, the round-3 numerics verifier's `v3_stage`; two runs each, identical; one case also at four threads) | 134 inputs: 67 returned EGH runs (767 features) and 50 returned Gaussian runs (445 features) beyond the earlier fixtures, isotope windows (`isowin` rows) across the averagine underflow and under NaN, tiny and full cutoffs, the charge wraps, the empty best pattern at `feature:min_isotope_fit` 0, `-0.0`, NaN (SIGSEGV) and `5e-324`, the `Precondition` order under ties, and the retention-time and intensity scales `vx_*`, `vy_*` (infinite widths and `FWHM` values: 7 of 9 at `1e37`, all from `1e38`) | `extended_cases_match_the_linux_release_build` |
| `fix4_stage` again (`../oracle/ffap-complete-fix5`, `run_onset.sh`, the round-4 numerics verifier's cases; two runs each, identical) | the onset of the width overflow: 20 retention-time scales and jittered inputs from `2e36` to `5e37`, Gaussian and EGH, plus the base case (first infinite width at `6e36`) | `width_onset_cases_match_the_linux_release_build` |
| `fix5_driver`, `ffap_instr_driver` (`../oracle/ffap-complete-fix5`; two runs each, identical but for the abort map's unique id) | the process-ending points outside the seed loop and a reused object's later terminations: 23 cases, exit status, the log's length on disk after each run and after the object is destroyed, every debug file's size and SHA-1 | `process_ending_refusals_outside_the_seed_loop_record_their_termination`, `a_reused_instance_leaves_the_executed_log_at_every_later_termination` |
| `fix4_reuse` (the verifier's `v3_reuse`; two runs each, identical) | one object run twice on FFC_1 in four scenarios: feature counts and the kept, appended and re-normalised isotope windows | `underflowed_windows_and_a_nan_cutoff_leave_nothing_to_append` |
| `fix4_vfi` (the instrumentation verifier's `vfi3_driver`; killed after 30 s, two runs each, identical) | a NaN retention time at scan 50 or 20 with `write_debug`, Gaussian and EGH, and without: the flushed `debug/log.txt` and the feature files of the seeds before the endless merge | `a_seed_loop_that_never_returns_keeps_what_the_executed_process_had_written` |
| the Release FeatureFinderCentroided on FFC_1 with the retention times scaled by `1e36` and `1e39` in the text (two runs each, identical apart from timing text) | exit 0 with `inf` intensities and, at `1e39`, `inf` widths in the featureXML | `infinite_feature_values_are_refused_by_the_featurexml_writer` (TOPP native difference 16) |

The FFC_1 score table also pins every loaded peak's m/z and intensity bits
against the C++ loader.

The oracle has two single-bin configurations (FFC_1 and the class test). The
B6 driver adds 7 and 10 bins, which exercise the four-cell interpolation, and up
to four charges; the degenerate captures add 1, 2 and 10 bins on degenerate
grids.

### The feature stage: what was compared and what it showed

`tests/feature_finder_picked.rs` reproduces each of the six executed
configurations twice: once through `run`, against the library's final
`FeatureMap` (tier 1), and once step by step through the public functions of
`extension.rs` and `fitting.rs`, against the C2 driver's per-seed replay
(adapted). 4,649 numeric values are compared in total.

**Every platform, since lead decision D10.** Since lane B3b the port's
Levenberg-Marquardt solver follows the Release build's Eigen kernels, and since
fix round 3 both fitters call the reference build's glibc 2.39 `exp` and `log`
as ported functions (`glibc_libm`: `__ieee754_exp_fma` and `__ieee754_log_fma`,
Arm optimized-routines with the executed FMA fusion, equal to the executed
library on 670 million probed inputs). Every compared value of the six
configurations is bit-identical to the Release capture on Linux x86_64 and on
macOS arm64: the five Gaussian configurations, seeds 11 and 12 of
`classtest_9247_tight_pattern` included (which departed by up to `2.3e-3`
from the product SDK), and the EGH configuration (`ffc1_asymmetric`), which
departed by up to `2.3038e-12` while its fit called the `libm` crate's `exp`,
`log` and `atan`. The numerics verifier of round 2 measured up to
`1.1096e-10` on further EGH configurations; those 16 configurations and 18
Gaussian ones are in `boundary_stage.tsv.gz` and are now exact too. The test
compares with `BITWISE` everywhere; the only other bound is `area_tolerance`,
for the `atan` of an EGH area on a host other than x86_64 Linux with the GNU C
Library (D10's fallback: the reference `atan` has no licence-clean upstream),
where the fixtures measured no departure on macOS arm64 (a measured maximum,
not a guarantee). Fix round 4 added 67 returned EGH runs (767 features) and
50 returned Gaussian runs (445 features) beyond these
(`extended_stage.tsv.gz`), all bit for bit on Linux x86_64 and macOS arm64.

Counts, charges, labels, `num_of_datapoints`, hull counts, hull point counts,
hull point coordinates, subordinate counts, abort reasons and abort counts are
compared exactly and agree everywhere. Every isotope-fit score, isotope-pattern
intensity and m/z score and every mass trace (peak identities, theoretical
intensities, baseline) is bit-identical on every platform.

**Platform notes (lead decision D8), what is left.** The fits no longer
depend on the host's `exp` and `log`; the earlier macOS arm64 bounds
(`5.4e-13` and the `1.1e-3` of seeds 11 and 12) and the unmeasured-platform
`1e-9` are gone. What remains host-dependent is:
- the `atan` of an EGH area on a host other than x86_64 Linux with the GNU C
  Library (the `libm` crate's), measured at 0 over these fixtures on macOS
  arm64, not a guarantee. The crate's `atan` differs from the reference's
  `__atan_fma` for 6.2 % of arguments in `[0, 10]`, 1.6 % with `|x|` in
  `[2^-14, 2^15)` and 0.02 % of random bit patterns (`2^28` inputs each,
  `../oracle/ffc-numerics-v3`); the `float` narrowing of the intensity hides
  nearly all of it (every fixture bit for bit on macOS arm64, the 67 returned
  EGH runs with 767 features of `extended_stage.tsv.gz` included), but the `double` area of `getArea` differs for
  a few percent of fits. A correctly rounded `atan` would depart about 90
  times less often (the reference misrounds 13 and 7 of 20,000 arguments of
  the first two ranges, the `libm` crate 1,219 and 324), but it is not the
  reference algorithm either and needs a new dependency; lead decision D13
  keeps the `libm` crate there, as a note for hosts other than the reference
  platform. On x86_64 Linux the host's `atan` is exact only where the GNU C
  Library selects `__atan_fma` (glibc 2.39 on a CPU with FMA, as on the
  reference node and the gate hosts); other glibc versions and CPUs are not
  measured;
- the sign of a NaN that the Levenberg-Marquardt iterations create from finite
  values (*Debug mode*, "The sign of a NaN inside a fit"), which no executed
  configuration reaches;
- the exactness of `f64::mul_add` on the host, which every Rust target
  provides as a correctly rounded fused multiply-add.

**The intended abundance override (adapted).** The library cannot compute
the override the source intends, so `../oracle/ffap-sem-completion/drivers/intended_abundance.cpp`
recomputes step 2.5 with a cleared override distribution, assigns the windows
to the protected `isotope_distributions_` and replays steps 3.1 to 4 with the
library's protected functions (two repetitions at one and four threads,
identical). For FFC_1 with `abundance_12C` 90 and 99 and `abundance_14N` 95 the
port's default reproduces every window bit for bit, the seeds with their
pattern and overall scores (18, 25 and 13), the candidates, the abort reasons
and the features (1, 8 and 2), bit for bit
(`the_intended_abundance_override_matches_the_adapted_release_replay`). The
executed library itself finds nothing in all three.

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
  - the refusals: charge wrap, abundances under
    `AbundanceOverride::Refuse`, and degenerate bin steps under
    `DegenerateBinStep::Refuse` only when the seed loop reads them; a NaN
    user-seed m/z among different ones and an unsorted spectrum with a NaN
    m/z, formerly refused, now give the executed outcome
    (`v3_seeds_nan_500_600`, `v3_mz_nan_three_unsorted`);
  - a single scan with a NaN intensity: sorted alone, scored with the zero
    steps' default NaN, and never reaching the seed loop (derived);
  - the conversions in `settings_follow_update_members`;
  - restriction and type violations;
  - a step-3.3.5 termination in a debug run keeps the seed's log lines and
    feature files (`a_step_3_3_5_termination_keeps_the_seed_debug_output`, a
    crate-internal unit test, since no executed input reaches it; fix round 5
    searched for one in vain, *Non-finite input*, refusal 6);
  - a window whose binary32 bins all underflow near `10^6` Da, a NaN trimming
    cutoff and a reused object's kept windows under it
    (`underflowed_windows_and_a_nan_cutoff_leave_nothing_to_append`; the
    underflow is derived from the exact bin values, far below the smallest
    binary32 subnormal);
  - the introsort on asymmetric and on other comparators
    (`asymmetric_comparisons_never_reach_the_guard`,
    `other_comparisons_are_refused_or_sorted_without_a_panic`).
- **Hand-derived (feature stage).**
  - `intersection_`: the two containment cases, both partial-overlap cases,
    disjoint and touching boxes, the division by the smaller total, and that a
    hull pair whose boxes do not intersect contributes nothing
    (`intersection_follows_the_source_cases`);
  - an isotope pattern that matched no peak is refused instead of dereferenced
    through the public function
    (`an_empty_pattern_is_refused_instead_of_dereferenced`; executed through
    `run`, see *Native differences*);
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
  - `score` at the grid centre, below the range (the wrapped `UInt` position)
    and at a NaN retention time (the NaN's sign cleared by the distance's
    absolute value).
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

The rows above predate combined fix round 3. Round 3 replaced the `libm`
crate's `exp` and `log` with the reference build's glibc 2.39 routines
(`glibc_libm`) and moved the fitter arithmetic to x86_64's operand order. It
was timed before (4ce747d) and after (e2a0661) with the same scratch test:
`run` at one thread, best of 9, `feature:rt_shape` symmetric (Gaussian
fits) and asymmetric (EGH fits), release build.

| Host | Symmetric, before / after | Asymmetric, before / after |
| --- | --- | --- |
| Linux, IBMI spock, AMD EPYC 9654 | 9.61 / 12.94 ms (+35 %) | 6.84 / 7.39 ms (+8 %) |
| macOS arm64 | 6.38 / 6.55 ms (+3 %) | 4.27 / 4.44 ms (+4 %) |

The Linux cost comes from the ported FMA fusion. The crate builds for the
baseline x86_64 target, which has no `fma` feature, so every `f64::mul_add`
in `glibc_libm` becomes an out-of-line `fma` call (Rust's
`compiler_builtins`, dispatched at run time).

Per-call timings over 2^24 inputs (`../oracle/ffap-complete-fix3/port-harness`):

| Host | Ported `exp` / `log` | glibc `exp` / `log` |
| --- | --- | --- |
| kim, baseline build | 26.9 / 27.6 ns | 4.5 / 4.0 ns |
| kim, `-C target-feature=+fma` | 3.3 / 3.3 ns | 4.5 / 4.0 ns |
| macOS arm64 (`mul_add` is one instruction) | 1.8 / 2.5 ns | 1.8 / 1.9 ns (Apple's) |

The `+fma` build matches every libm probe. Neither of the two ways to remove
the cost fits in this lane:

- Building with `-C target-feature=+fma` is a build-configuration decision.
  It also leaves out CPUs without FMA.
- Dispatching at run time needs `unsafe`, which the crate forbids.

**Test time (fix round 4).** Two cases of `extended_stage.tsv.gz`
(`vw_w1`, `vw_w1_ipo0`) compute 274,001 isotope windows each, twice (the
source's run and the `DegenerateBinStep::Refuse` rerun). An unoptimised test
build takes about 200 µs per heavy window on macOS arm64, so each of those
cases takes close to two minutes; `extended_cases_match_the_linux_release_build`
therefore replays its cases on up to eight threads and takes 129 s there
instead of 276 s serially. The release build computes the same windows in a
fraction of that. Lead decision D13 accepts the longer
`feature_finder_picked` target (about three minutes on the gate hosts); no
assertion is dropped to shorten it. Fix round 5 adds 21 short stage cases
(`width_onset_cases_match_the_linux_release_build`) and two instrumentation
tests of 23 executed cases, each a few seconds, run on up to four threads.

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
  --test isotopes --test isotopes_source_precision --test mass_trace --test mass_trace_detection \
  --test topp_feature_finder_centroided --test feature_finder_picked_instrumentation \
  --test progress_logger --test topp_threads --test parallel_determinism \
  --test param --test default_param_handler
cargo +1.85.0 check --locked --all-features --all-targets
cargo test --locked --all-features --lib feature_finder_picked
cargo test --locked --all-features --lib isotopes
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
6. **`FeatureFinderDefs` is defined in two headers.**
   `FEATUREFINDER/FeatureFinderAlgorithmPicked.h:24-55` and
   `FEATUREFINDER/FeatureFinderDefs.h:19-50` define the same `struct
   OPENMS_DLLAPI FeatureFinderDefs`. A translation unit that includes both does
   not compile: GCC 14 reports `redefinition of 'struct
   OpenMS::FeatureFinderDefs'` (executed with the Release install's compiler,
   `../oracle/ffap-sem-completion/drivers/defs_both_headers.cpp`). Nothing
   includes `FeatureFinderDefs.h` today, so the duplicate is latent.
7. **Short inputs pass an inverted range to `startProgress`.** For an input
   with fewer than `2 * min_spectra_` scans, steps 2 and 3.2 call
   `startProgress(min_spectra_, n - min(min_spectra_, n))` (`.cpp:297-298`,
   `:493-494`) with `begin > end`. The Debug build stops there with
   `OPENMS_PRECONDITION(begin <= end)` (`ProgressLogger.cpp:235`, exit 8 for
   `FileFilter_44_input.mzML` with `-force` and for
   `FileConverter_31_output.mzML`); the Release build runs on and returns an
   empty map, because the seed loop is empty for every input of at most
   `2 * min_spectra_` scans. Neither build tells the user that the input is too
   short to hold a seed. This is the executed content of the zero-width record
   `CPP-312`, whose empty map came from the short input, not from the zero
   extent.

8. **Step 1 drops scans with a non-finite drift time without a word.**
   `areaBeginConst` sets the mobility range to `RangeMobility{}` made
   non-empty, `[lowest, max]` (`MSExperiment.cpp:562-571`), and
   `AreaIterator::nextScan_` skips every scan outside it
   (`AreaIterator.h:277-298`). A scan whose drift time is NaN or infinite
   contributes no intensity to the quantiles, although its peaks are scored
   against them and can become seeds (executed: `v2_dt_*`, `v3_dt_*`; every
   drift time NaN leaves every cell empty and finds 10 features instead of 8).
   Proposed fix: iterate without a mobility filter here, or refuse such input.
9. **`feature:min_isotope_fit` 0 crashes the seed loop.** `findBestIsotopeFit_`
   returns 0 and leaves `best_pattern` empty when no placement qualifies, and
   `isotope_fit_quality < min_isotope_fit_` is false at 0, so
   `extendMassTraces_` reads `pattern.spectrum[0]` of an empty vector
   (`.cpp:1347-1349`). Executed: SIGSEGV on FeatureFinderCentroided_1 with
   `feature:min_isotope_fit 0` and `seed:min_score 0` (library and tool, twice
   each), and on three negative-intensity inputs. Proposed fix: skip a seed
   whose pattern is empty, or make the parameter's minimum positive.
10. **Heavy averagine windows are silently empty.** From an averagine mass
    of about 273,770 Da on (executed between 273,769.5 and 273,770.5 Da; the
    windows just below keep a single bin), all 20 binary32 bins of
    `estimateFromPeptideWeight` underflow;
    `IsotopeDistribution::renormalize` divides by the zero sum and `trimRight`
    discards the NaN weights, so every peak of such a mass gets pattern
    score 0 and no feature, without a message (executed: m/z `136,850.5` to
    `1e6` at charge 2). Proposed fix: compute the averagine in `double` or in
    log space, or report the limit.
11. **A NaN `isotopic_pattern:intensity_percentage_optional` is accepted.**
    `ParamEntry::isValid`'s range comparisons are false for NaN, and the NaN
    cutoff then empties every window: no feature, no message (executed:
    `p_ipo_nan*`). The same holds for the other percentage and tolerance
    parameters the verifiers executed. Proposed fix: reject NaN in
    `isValid`.
12. **Unsorted input with a mis-sized data array half-sorts the caller's
    map.** `run(PeakMap&&, ...)` sorts the caller's object in place, and
    `sortSpectra` throws `Exception::Precondition` from `MSSpectrum::sort`
    after it has reordered the spectra and sorted the earlier ones' peaks, so
    the caller's map is left partly sorted (executed: the `a_*` cases of
    `boundary_stage.tsv.gz`). Proposed fix: validate the arrays before
    sorting anything.
13. **The step-1 progress range wraps.** `startProgress(0, intensity_bins_ *
    intensity_bins_)` multiplies two `UInt`s, so 65,536 bins announce a range
    of 0 (executed); `setProgress` then reports values beyond the end.
    Cosmetic. Proposed fix: multiply in `Size`.
14. **Large retention times silently give infinite feature widths.**
    `setWidth(fitter->getFWHM())` narrows the `double` FWHM to `float`
    (`.cpp:742`), and the intensity `getArea() / max` too (`.cpp:790`). With
    FFC_1's retention times scaled and its parameters unchanged, the first
    infinite width appears at a scale of `6e36` (1 of 9 features; 3 of 9 at
    `8e36`, 7 of 9 at `1e37`, all from `1.5e37`; none at `4e36`), and every
    intensity is infinite from `1e36` on (finite at `1e33`): the run returns
    features with an infinite width, `FWHM` meta value and intensity and
    writes `inf` into the featureXML, without a message (executed: `nb_rt_*`,
    `vy_rt_1e36` to `vy_rt_1e150`, and the Release FeatureFinderCentroided on
    FFC_1's input with every scan start time times `1e36` and `1e39`).
    Proposed fix: check the fitted width and area against the `float` range,
    or store them as `double`.

### From the seed stage (B6)

Items 1 and 2 are executed; the others come from source review.

1. **`mass_trace:min_spectra = 1` silently finds nothing.** `min_spectra_ =
   floor(1 * 0.5) = 0`. Every trace score becomes 0/0 (`.cpp:340`), every peak a
   local maximum, and every overall score NaN. The run reports 0 seeds and
   0 features without a warning (executed: B6 driver `ffc1_min_spectra_1`).
   The parameter's minimum should be 2, or the value should be rejected.
2. **The overall score depends on the platform `powf`.**
   `std::pow(float, float)` at `.cpp:506` calls the C library `powf`, which is
   not correctly rounded: 99 of 30,840 executed scores are one binary32 step off
   with Apple's (macOS arm64 product SDK) and 8 with glibc 2.39's (Linux x86_64
   Release build, FMA variant), different scores in each case, so a score next
   to `seed:min_score` can flip a seed between platforms. Evaluating in
   `double` and rounding once gives the correctly rounded value. (The port
   reproduces the reference build's glibc `powf`, lead decision D5.)
3. **`write_debug` throws.** `writeFeatureDebugInfo_` reads
   `debug:pseudo_rt_shift` (`.cpp:2137`), but the declared parameter is
   `advanced:pseudo_rt_shift` (`.cpp:124`). The resulting `ElementNotFound`
   escapes the OpenMP region. Executed: the Release tool dies of SIGABRT (shell
   status 134) in every such run (tool cases a4, b1, c1).
4. **Degenerate intensity bins are undefined behaviour.** With one RT or one
   m/z, or an extent that underflows in the division, `intensity_rt_step_` or
   `intensity_mz_step_` is 0, and an overflowing RT extent makes it infinite;
   `intensityScore_` then converts `floor(NaN)` or `floor(inf)` to `UInt`
   (`.cpp:1837-1838`). The Linux x86_64 Release build (executed, see
   *Degenerate intensity bins*) makes every intensity score NaN and finds
   nothing, silently: an input whose retention times or m/z values are all
   equal yields an empty map with no message. The steps should be checked and
   the input refused, or the intensity score defined for a single bin.
5. **The charge count wraps.** `UInt charge_count = charge_high -
   charge_low + 1` (`.cpp:197`) wraps for `charge_low > charge_high + 1`, and
   `3 + 2 * charge_count` wraps too, also for `charge_low` 1 and
   `charge_high` `INT_MAX`; the float arrays are then written past their end,
   or allocated at nearly 2^32 per spectrum (executed: SIGSEGV for 4/2,
   `INT_MAX`/1, `INT_MAX`/498 and 1/`INT_MAX`; `std::bad_alloc` for 5/2, 6/2,
   7/2 and 2/`INT_MAX`). Proposed fix: check `charge_low <= charge_high` and
   bound the count.
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

- **What is ported.** Every declaration of `FeatureFinderAlgorithmPicked.h` has a
  counterpart or a recorded reason in the API mapping above. That includes
  `FeatureFinderDefs` (`defs.rs`) and the whole algorithm, with the degenerate
  intensity bins reproduced as the Linux x86_64 Release build computes them
  (port/ffap-semantics), and `writeFeatureDebugInfo_`, `abort_reasons_`, the
  `ProgressLogger` base and the whole `write_debug` mode, together with the
  instance state and `run` into a caller's map (port/ffap-instrumentation).
  What stays refused there is where the source terminates or is not
  reproducible (*Debug mode*, *Reusing an instance*).
- **What is refused, and why each refusal is as narrow as its reason.**
  - `DegenerateBinStep::Refuse` is an opt-out, not the default; it refuses only
    a zero or infinite step whose scores the seed loop reads.
  - The `UInt` score-array count `3 + 2 * charge_count` that wraps
    (`charge_low > charge_high + 1`, and `charge_low` 1 with `charge_high`
    `INT_MAX`), refused where the count is computed, whatever the `Limits`
    (lead decision D12): an out-of-bounds write for the counts `2^31 - 1`,
    `-1` and `-4` and below (executed SIGSEGV, recorded as a `ScoreArrays`
    termination), and for `-2` and `-3` an in-bounds allocation of about 2^32
    arrays per spectrum whose outcome depends on memory (executed
    `std::bad_alloc`), refused unconditionally by lead decision D13. The
    in-bounds wrap of the step-1 progress range is reproduced.
  - Non-finite input is read as the Release build reads it (189, 52, 85 and
    134 executed cases). What remains refused is listed in *Non-finite input*: the
    endless NaN profile merge, the out-of-bounds read of an empty best
    pattern (reachable with `feature:min_isotope_fit` 0, executed SIGSEGV),
    the wrapped charge count and the step-3.3.5 exception that terminates the
    source; the introsort's out-of-bounds guard is kept but unreachable for
    the algorithm's comparisons. Each of them records its `DebugTermination`
    (the seed-loop ones after the seed's log lines, the charge count before
    the log opens), with or without `write_debug`, with the length at which
    the executed process leaves `debug/log.txt`.
  - On a caller's map or a reused object: a charge-0 feature in an
    overlapping pair of different charges (step 4's `%` traps, executed
    SIGFPE) and a stale abort seed outside the current input (executed
    SIGSEGV), both recorded as terminations (*Reusing an instance*).
  - Averagine windows whose binary32 bins all underflow and a NaN
    `intensity_percentage_optional` are reproduced (empty windows), not
    refused.
  - 64-bit integer parameters are narrowed and converted as the Release build
    does (executed), not refused; only the source's own `InvalidParameter` and
    `ConversionError` remain, with their texts.
  - The `Limits` ceilings (bounded work); step 2.5 above `vector::max_size()`
    returns the source's `length_error` text instead.
  - Non-finite and negative FWHM, `score_fit`, `score_correlation` and
    `EGH_*` values are stored as the source stores them since fix round 4
    (an infinite width is reachable from retention times of `6e36` and more,
    executed), not refused, through the crate-private
    `MetaValue::source_float` that lead decision D13 accepts; the crate's
    checked consumers (`validate`, the featureXML writer) still refuse such
    features. The featureXML writer's refusal and the CLI's "Unable to read
    file" wording for that write failure are split off into a separate task
    (D13); TOPP native difference 16 stays as recorded.
  - `AbundanceOverride::Refuse` is an opt-out; the default computes the
    intended override, the one designed difference (`CPP-247`), now pinned
    against an adapted Release replay.
- **Unreachable source behaviour.** `Exception::UnableToFit` cannot be thrown
  from the seed loop (argument at `FittedModel::fit`, including the `int`
  residual count), and the `aborts_` data race at more than one thread has no
  reproducible C++ result, so neither is tested against the C++. The race on
  `aborts_`, `abort_reasons_` and `log_` is answered with the single-thread
  result, not refused: lead decision D11 accepts it as the one documented
  exception to D1, because the determinism contract requires parallel output
  to equal serial output and a refusal would block essentially every parallel
  run.
- **The C library.** Every libm transcendental on the path is the reference
  build's (lead decision D10): `powf` (D5), `exp` and `log` ported from Arm
  optimized-routines with the executed FMA fusion, bit for bit on every
  platform; `atan` has no licence-clean upstream of its algorithm (IBM's, LGPL
  only), so it is the host's with glibc (exact on the reference platform) and
  the `libm` crate's elsewhere, where the EGH area's bound is a measured
  maximum; `sqrt` is correctly rounded everywhere.
- **Lead decisions D13 (after fix round 4).** The crate-private
  `MetaValue::source_float` is accepted (the public API still refuses
  non-finite values); the featureXML writer's refusal of non-finite values
  and the CLI's wording for a write failure are a separate task, and TOPP
  native difference 16 stays; the longer `feature_finder_picked` test time is
  accepted; the charge counts `-2` and `-3` stay refused unconditionally; the
  `libm` crate's `atan` off x86_64 Linux with glibc stays, as a note for
  non-reference platforms; the step-3.3.5 termination keeps its
  source-review status after one more search for an executed input
  (*Non-finite input*, refusal 6).
- **Reproduced undefined and unspecified behaviour** (lead decisions D1 to D3):
  the float-to-integer conversions, the binary searches and both sorts on NaN
  keys, the order of equal sort keys, and the address-independent text of a
  string or list shift, each where its outcome is measured, repeatable and
  explained by the executed instructions and every read stays in bounds.
- Rust files: `src/analysis/feature_finder_picked/algorithm.rs`, `scoring.rs`,
  `seeds.rs`, `extension.rs`, `fitting.rs`, `resolution.rs`, `defs.rs`,
  `instance.rs`, `debug.rs`, `source_sort.rs`, `glibc_powf.rs` and
  `glibc_libm.rs`. Tests:
  `tests/feature_finder_picked_seeds.rs`, `tests/feature_finder_picked.rs`,
  `tests/feature_finder_picked_instrumentation.rs` and
  `tests/topp_feature_finder_centroided.rs`.
- The ledger scope's "steps 3.3 to 5" should read "steps 3.3 and 4": the
  source's last step is step 4 (`.cpp:860`).
- `docs/doc-coverage.json` needs `--write` (every module of the directory is at
  100 %; `glibc_powf.rs` and `glibc_libm.rs` are crate-private); the floor was
  not re-recorded here because the file is the integrator's.
