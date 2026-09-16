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
| [`algorithm.rs`](../src/analysis/feature_finder_picked/algorithm.rs) | defaults, `Settings`, `validate_input`, `Limits`, `Options` with `AbundanceOverride` and `DegenerateBinStep`, `PseudoRtShiftKey`, `RejectedParameters`, `run` and `run_with_options` (a fresh map), `feature_stage` (stage level) |
| [`scoring.rs`](../src/analysis/feature_finder_picked/scoring.rs) | `position_score`, `nearest_from`, `IntensityThresholds`, `ScoreArrays`, `find_isotope`, `isotope_score`; the crate-private `x86_64` emulation of the Release build's `UInt` conversion and NaN rules |
| [`defs.rs`](../src/analysis/feature_finder_picked/defs.rs) | `FeatureFinderDefs`: `IndexPair`, `IndexSet`, `ChargedIndexSet`, `Flag`, `NoSuccessor` |
| [`seeds.rs`](../src/analysis/feature_finder_picked/seeds.rs) | `IsotopeWindows`, `SeedStage`, `ChargeSeeds`, `UserSeed`, `overall_score` |
| [`extension.rs`](../src/analysis/feature_finder_picked/extension.rs) | `OverallScores`, `find_best_isotope_fit`, `extend_mass_traces`, `extend_mass_trace` |
| [`fitting.rs`](../src/analysis/feature_finder_picked/fitting.rs) | `FittedModel`, `crop_feature`, `check_feature_quality`, `build_feature`, the abort reasons |
| [`resolution.rs`](../src/analysis/feature_finder_picked/resolution.rs) | `intersection`, `resolve_overlaps`, `annotate_apex` |
| [`instance.rs`](../src/analysis/feature_finder_picked/instance.rs) | `FeatureFinderAlgorithmPicked`, the source object with its state across runs (parameters, seeds, aborts, abort reasons, the log stream, the isotope windows, progress), and its `run` into the caller's map |
| [`debug.rs`](../src/analysis/feature_finder_picked/debug.rs) | `DebugOutput`, `DebugLog`, `SeedMap`, `AbortReasons`, `seed_map`, `abort_map`, `debug_experiment`, `write_feature_debug_info`, `PseudoRtShift`, `HEAP_ADDRESS_END`, `ReportLine` |
| [`source_sort.rs`](../src/analysis/feature_finder_picked/source_sort.rs) | `source_sort_by`, `source_sort_reversed_by`, `source_sort_permutation`: the order in which the Linux x86_64 Release build's `std::sort` (libstdc++ introsort) leaves equal elements |

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
(steps 0 to 3.2, the degenerate intensity bins and `FeatureFinderDefs`),
[`tests/feature_finder_picked.rs`](../tests/feature_finder_picked.rs)
(steps 3.3 and 4, and the intended abundance override) and, through the tool,
[`tests/topp_feature_finder_centroided.rs`](../tests/topp_feature_finder_centroided.rs);
[`tests/feature_finder_picked_instrumentation.rs`](../tests/feature_finder_picked_instrumentation.rs)
(the instance, debug mode, progress and a caller's map; manifest
[`tests/data/feature_finder_picked_instrumentation_provenance.json`](../tests/data/feature_finder_picked_instrumentation_provenance.json),
which needs `featurexml` too).
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
| base `DefaultParamHandler` | `instance::FeatureFinderAlgorithmPicked`: `name`, `set_name`, `subsections`, `handler_equal` (`operator==`), `defaults`, `default_parameters`, `parameters`, `set_parameters`, `set_parameters_logged`, over the crate's `DefaultParamHandler`; `Settings::from_parameters` for the stateless `run`. The static `writeParametersToMetaValues` is `DefaultParamHandler::write_parameters_to_meta_values`; the protected `check_defaults_` and `warn_empty_defaults_` are not settable, as the source object never changes them | `setName` renames the handler in the unknown-parameter warnings, `getSubsections` is empty, and `operator==` compares only the handler part, with the refused set the source keeps in `param_` (source-reviewed, `the_handler_base_renames_compares_and_has_no_subsections`). Executed: `getParameters` after construction, after a run and after an empty-input run; the unknown-parameter warnings in `checkDefaults` order, and on a refused set only those before the refused entry; after a refused set `parameters()` shows that set merged with the defaults while the settings stay, as the source assigns before it checks (`algorithm::RejectedParameters::Shown`, the default; `Discarded` keeps the accepted set) |
| base `ProgressLogger` | `instance::FeatureFinderAlgorithmPicked::set_log_type`, `log_type`, `set_progress_logger`, `progress_logger_mut`; the base's public `startProgress`, `setProgress`, `nextProgress` and `endProgress` are the methods of the logger `progress_logger_mut` returns, and are absent under type `NONE` (where the source's calls do nothing), as `progress_logger_mut` then returns `None` | the 20 call sites of `.cpp:241-992`: every start, set and end call with its label, range and value, equal to the Release build's complete event sequence (executed with a counting clock, see *Debug mode*); no logger (type `NONE`, the default) costs nothing. An inverted range (steps 2 and 3.2 on fewer than `2 * min_spectra` scans) is passed with `end = begin`, the one recorded difference; the `CMD` output is the same |
| `MapType`, `SpectrumType`, `FloatDataArrays` | `MSExperiment`, `MSSpectrum`, `ScoreArrays` | see *Score arrays* below |
| `PeakType`, `Seed`, `MassTrace`, `MassTraces`, `TheoreticalIsotopePattern`, `IsotopePattern` (protected) | `Peak1D` and the `helper_structs` types | |
| `FeatureFinderAlgorithmPicked()` | `instance::FeatureFinderAlgorithmPicked::new`, `with_options`; `default_parameters`, `HANDLER_NAME` | |
| `setSeeds(const FeatureMap&)` | `instance::FeatureFinderAlgorithmPicked::set_seeds`, `seeds`; `run` replaces them, as the source's `run` does; the `seeds` argument of `SeedStage::run` | |
| `setData_(MSExperiment&&, FeatureMap&)` (private) | `instance::FeatureFinderAlgorithmPicked::run` consumes the experiment and extends the caller's `&mut FeatureMap` | |
| `run(PeakMap&&, FeatureMap&, const Param&, const FeatureMap&)` | `instance::FeatureFinderAlgorithmPicked::run`, which extends the caller's map as the source does (see *Reusing an instance*); `algorithm::run` and `run_with_options` for a fresh map | the stateless functions return `RunOutput` with the features, the log, the abort counts and the debug output |
| `getDefaultParameters()` | `default_parameters()` | 29 entries, 7 section descriptions |
| `run_()` (protected) | `instance` `run_core`: `SeedStage::prepare` (steps 0 to 2.5), then per charge `select_next_charge`, `extend_charge` and `settle_charge`, then step 4; `SeedStage::compute` and `feature_stage` at stage level | |
| `map_` | `SeedStage::experiment`, `into_experiment` | |
| `features_` | the caller's `&mut FeatureMap` of `instance::FeatureFinderAlgorithmPicked::run`; `RunOutput::features` of the stateless `run` | |
| `log_`, `debug_` | `Settings::write_debug`; `debug::DebugOutput` (`log`, `log_opened`) returned by `run` | the log text in source order and formatting, with the stream's state across runs (see *Debug mode*); the tool writes the file |
| `aborts_`, `abort_()` | `instance::FeatureFinderAlgorithmPicked::aborts` (`u32`, accumulating over runs); `RunOutput::aborts` | counted serially in seed order: the source's single-thread counts at any thread count; the source races them (candidate 5) |
| `abort_reasons_` | `instance::FeatureFinderAlgorithmPicked::abort_reasons`, `debug::AbortReasons` (keyed by intensity, as `Seed::operator<`), `debug::abort_map` | never cleared, as in the source; the abort map reads the stored indices in the current input (see *Reusing an instance*) |
| `seeds_` | `SeedStage::user_seeds` (`UserSeed`: m/z and RT) | the only fields the source reads |
| `pattern_tolerance_` ... `max_feature_intersection_`, `reported_mz_` | `Settings` fields of the same names; `reported_mz` as `ReportedMz` | |
| `intensity_rt_step_`, `intensity_mz_step_`, `intensity_thresholds_` | `IntensityThresholds::rt_step`, `mz_step`, `quantiles` | |
| `isotope_distributions_` | `IsotopeWindows::patterns`; `IsotopeWindows::precalculate_onto` and `instance::FeatureFinderAlgorithmPicked::isotope_windows` keep them across runs | extended, never cleared, as in the source (see *Reusing an instance*) |
| `updateMembers_()` | `Settings::from_parameters` | also reads the `run_` locals: charges, `fit:max_iterations`, abundances, seed thresholds, user-seed tolerances, `feature:min_score`, `feature:rt_shape` |
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
| `chooseTraceFitter_(double&)` | `fitting::FittedModel::new` from `Settings::rt_shape` | the enum replaces the pointer plus the `tau != 0` flag and the `dynamic_pointer_cast` |
| `cropFeature_()` | `fitting::crop_feature` | returns the cropped traces instead of an out-parameter |
| `checkFeatureQuality_()` | `fitting::check_feature_quality` | returns `QualityOutcome`: `Accepted(FeatureQuality)` or `Rejected(reason)` |
| step 3.3 (`.cpp:576-856`), the `omp parallel for` and the containment pass | `instance::extend_charge` and `settle_charge`, with `fitting::build_feature` for step 3.3.5 | `concept::parallel::map_collect` over the seed indices |
| step 4 (`.cpp:859-1016`) | `source_sort::source_sort_by` (the Release build's introsort order), `resolution::resolve_overlaps`, `retain`, `source_sort_by` with the arguments swapped, `resolution::annotate_apex` | |
| `writeFeatureDebugInfo_()` | `debug::write_feature_debug_info`, `debug::FeatureDebugFiles`, `debug::PseudoRtShift`; `algorithm::PseudoRtShiftKey` | the `.dta`, `_cropped.dta` and `.plot` texts, byte-identical to the Release build; the undeclared `debug:pseudo_rt_shift` is read as the source reads it, a string or list value included (see *Debug mode*) |
| `operator=`, copy constructor (private, not implemented) | not applicable | the stage is an owned value |
| step 1, second half (`.cpp:280-287`) | `fill_intensity_scores` (crate-private) | |
| step 2 (`.cpp:291-348`) | `fill_trace_scores` (crate-private) | |
| step 3.1 (`.cpp:447-488`) | `fill_pattern_scores` (private) | |
| step 3.2 (`.cpp:489-574`) | `select_seeds` (private), `overall_score`, `SeedStage::charges` | |

Native additions: `Limits`, `Options`, `AbundanceOverride`, `DegenerateBinStep`, `RunOutput`,
`validate_input`, `UNSORTED_WARNING`, `BASE_MAX_ISOTOPES`,
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
`SIZE_MAX`).

## Debug mode

`write_debug = true` makes the source write into the working directory
(`FeatureFinderAlgorithmPicked.cpp:226-232`, `550-571`, `714`, `1028-1052`).
The library port returns the same content as data in `debug::DebugOutput`
(`RunOutput::debug`, `FeatureFinderAlgorithmPicked::debug_output`), and
`FeatureFinderCentroided` creates `debug/` and `debug/features/` and writes the
files under the source's names:

| Source output | `DebugOutput` | Evidence (Linux x86_64 Release, `OMP_NUM_THREADS=1`) |
| --- | --- | --- |
| `debug/log.txt` (`log_`, 58 write statements) | `log` (`DebugLog`), `log_opened` | byte-identical for the tool cases a1, a2, a3 and the driver's two-run object; `double` values printed as `operator<<` prints them, glibc's `nan`/`-nan` included |
| `debug/seeds_<charge>.featureXML`, per charge, also for a charge without seeds | `seed_maps` (`SeedMap`) | D6 (decoded, ids excluded): a1, a2, a3, a4, b1, stale scaled, debug_twice |
| `debug/features/<plot_nr>.dta`, `_cropped.dta`, `.plot` (`writeFeatureDebugInfo_`) | `feature_files` (`FeatureDebugFiles`, `debug::write_feature_debug_info`) | byte-identical: the 75 files and the log of each of the driver cases declared-shift500, -shift123, -int250, -egh and -prefilled, of a string and a string-list shift, and of a string shift with a scan at RT `5e-275` and `1e-289` (`debug_digests.tsv`) |
| `debug/abort_reasons.featureXML` | `abort_reasons` (`debug::abort_map`) | D6: a1, a2, declared-*, debug_twice, stale scaled; the feature ids `0, 1, ...` exactly |
| `debug/input.mzML`: the input with the score arrays, without the overall score | `input` (`debug::debug_experiment`) | D6: a1, a2, a3, debug_twice; NaN scores compared as NaN. `mzml::write_source_float_arrays` writes the non-finite values the source writes |
| the process terminates in `writeFeatureDebugInfo_` | `termination` (`DebugTermination`) and `Error::Unsupported` | a4 and b1: charge, exception and message; `log.flushed_bytes()` is the length of the executed file after the SIGABRT |

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
  retention time of `2^-974` (about `6.3e-294`; `2^-973` for trace 2) and at
  a fitted centre between about `1.4e-304` and `1.4e-303`, depending on its
  digits. A bound from the whole 56-bit address space would refuse `1e-289`,
  which the Release build writes reproducibly.
- `PseudoRtShiftKey::Declared`: read `advanced:pseudo_rt_shift` and write the
  feature files for every seed that reaches the fit, which is what the source
  evidently intends. The driver cases declared-* (with `debug:pseudo_rt_shift`
  set, so the Release build completes) give the same bytes under both
  policies.

A fit that throws `UnableToFit` (`TraceFitter.cpp:111`, `:129`) inside the
same parallel loop would end the source process at `.cpp:670`, before the
debug write at `:714`, but no input reaches either throw (see the
`Exception::UnableToFit` row of *Native differences*). An error from the
port's own fit (one of its ceilings, or the refused NaN profile merge) ends
the run with that error at the first such seed in seed order; no executed
configuration reaches a failing fit.

**The file at termination.** `log_` is an `std::ofstream` whose 8,191-byte
buffer is lost when the process aborts. `DebugLog` models libstdc++'s
`basic_filebuf` (`fstream.tcc`: a block write when the text does not fit in
the free space, one `sputc` per `char` insertion) and records how many bytes
reached the file. The model predicts the executed files of a4, b1 and of the
two-run object while it was still alive; the tool writes that prefix.

**Order.** The source writes the log from inside its parallel loop, under a
critical section, in schedule order. The port collects every seed's lines and
aborts and appends them in seed order: the source's single-thread output at
every thread count (`debug_output_is_identical_at_every_thread_count`). The
executed tool at four threads wrote a different log in each of three
repetitions (cases c1, c2); that output is undefined and not compared.

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
(`the_progress_event_sequence_matches_the_release_build`). The one difference
is the inverted range the source passes in steps 2 and 3.2 on an input of
fewer than `2 * min_spectra` scans (`S 5 0`): the port's
`ProgressLogger::start_progress` refuses an inverted range, so the port passes
`end = begin` (`S 5 5`). No `setProgress` call falls into such a step, so the
`CMD` output is the same.

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
- **`log_` is opened once.** A second debug run's `open` fails on the open
  stream, the stream's `failbit` drops every write, and the file keeps the first
  run's text (`DebugOutput::log_opened` is `false`).
- **Parameters.** `getParameters()` is the last run's merged set, or the
  last refused set: `setParameters`, which `run` calls first, assigns the new
  set before it checks it, and a refused set stays visible while the members
  keep the accepted values (`rejected_stdout.txt`: a direct call and a run).
  The warnings `checkDefaults` logs before it throws are the unknown keys it
  visits before the refused entry; a run puts them in its report.

Two undefined continuations are refused at the point where the source's
behaviour ends: a stale abort seed whose indices lie outside the current input
(`Error::InvalidValue` while building the abort map; the executed driver died
of SIGSEGV in all seven repetitions) and a caller's feature of charge 0 in an
overlapping pair with a different charge (`Error::InvalidValue`; the source's
`%` traps, SIGFPE in both repetitions). A NaN key that made the source's
`std::sort` read outside the map would be refused as well; no executed input
did (`sort_keys.txt`, three NaN sequences).

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

The source reads infinite and NaN retention times, m/z values, intensities and
user-seed positions without a check. The port follows it, and refuses only
where the source's behaviour has no reproducible answer. What the executed
Linux x86_64 Release build does was measured with the driver
`nonfinite_stage` (`../oracle/ffap-sem-completion/drivers/nonfinite_stage.cpp`)
on 189 modified FeatureFinderCentroided_1 inputs, each run twice (identical;
five also twice at four threads, identical apart from the one-thread abort
rows). The fixture `nonfinite_stage.tsv.gz` holds every outcome, and
`non_finite_inputs_match_the_linux_release_build` replays all 189.

**What the source does, and what the port does.**

- *The input check* (`run`, `.cpp:1057-1098`). `MSExperiment::isSorted(true)`
  compares retention times with `>` and each spectrum with `std::is_sorted`;
  both are false for a NaN, so a NaN never makes the input unsorted. `-inf`
  sorts before every m/z and fails the positive-m/z check (`IllegalArgument`).
  The port runs the same comparisons (`validate_input`) and sorts as
  `sortSpectra` does.
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
  the same cell contents. An infinite intensity sorts to the end of its cell
  and changes the quantiles; the score of that peak is `inf / inf`, NaN.
- *Step 2.* The trace scores compare intensities only, so an infinite
  intensity changes local maxima; a NaN or infinite m/z difference fails both
  comparisons of `positionScore_` (score 0). `findNearest` is reproduced with
  the same `lower_bound`.
- *Step 2.5.* `Size num_isotopes = ceil(max_mass / width) + 1` is converted by
  `comisd`/`cvttsd2si`/`btc` (`libOpenMS.so` `0x18e46f4`-`0x18e46fe`,
  `0x18e6c3b`-`0x18e6c44`): an infinite maximum m/z, or a count of `2^64` or
  more (m/z `1e300`), gives **no** window. A count in `[2^63, 2^64)` makes
  `resize` throw `std::length_error` (`vector::_M_default_append`); the port's
  `Limits::max_isotope_windows` refuses it first. The port converts with the
  crate-private `x86_64::truncate_to_u64`, which reproduces the sequence.
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

**Measured outcomes** (FFC_1 with its INI; `in`/`innear` sweeps set one peak of
or next to each of the 25 seeds):

| Cases | Executed Release outcome | Port |
| --- | --- | --- |
| `+inf`/`-inf` intensity at 3 fixed positions and 136 sweep positions (also with `seed:min_score` 0, `mass_trace:min_spectra` 2, the EGH model, `feature:min_isotope_fit` 0, reported m/z `average` and `maximum`) | 7 to 13 features; every score, seed, feature, meta value and hull recorded | the same (bit for bit on Linux x86_64; see below for macOS) |
| `-inf` m/z (first peak, last peak, a whole spectrum) | `IllegalArgument`, positive m/z | the same text |
| `+inf` m/z (last peak, a whole first spectrum, 5 spectra), m/z `1e300` | `InvalidValue`, no window | the same text |
| m/z `6.9e20` (`2^63 <= count < 2^64`) | `std::length_error` | `Limits::max_isotope_windows` |
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
| NaN intensity (2 fixed positions, 6 sweep positions); NaN seed m/z among others; NaN retention time with an unsorted scan | 8, 7, 7 features; never returns | refused at the sort (below) |

**What is refused, and where.** Each refusal is `Error::InvalidValue`, at the
first point where the source becomes undefined:

1. *A NaN key in a `std::sort` or `std::stable_sort`* whose other keys are not
   all bit-identical to it: the bin intensities of step 1 (`.cpp:270`), the
   spectra by retention time and the peaks of an unsorted spectrum by m/z
   (`MSExperiment::sortSpectra`), the chromatograms by product m/z and the
   peaks of an unsorted chromatogram by retention time
   (`sortChromatograms`), and the user seeds by m/z when two of their other
   m/z values differ (`.cpp:190`). The standard requires a strict weak
   ordering there, which a NaN breaks as soon as two other keys differ; when
   every other key is equivalent, the order of the equivalent elements is
   unspecified and observable (it decides the stored quantiles or the spectrum
   order). The Release build's order is libstdc++'s introsort order, which
   `source_sort` reproduces for the seeds and the feature map but which these
   sites do not use yet. Executed: `int_nan_50_3` (8 features), `int_nan_40_10` (7),
   `seeds_mz_nan` (7) and the six `sw_nan_next2_*` returned; the fixture keeps
   their outcomes. The user seeds are exempt when every other m/z is one value,
   because `near_user_seed` then answers the same for every order.
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
4. *A feature m/z without an isotope window at step 3.3.5* (`.cpp:790`): the
   source's `InvalidValue` leaves its OpenMP region uncaught and
   `std::terminate` ends the process. A NaN feature m/z (an infinite intensity
   kept in the reported traces) would reach it; no executed input did
   (source review).

**The platform split.** On Linux x86_64 with glibc every returned feature is
bit for bit. On macOS arm64 the Gaussian fits split as in
`docs/TRACE_FITTER_SUPPORT.md`: of the 1,135 compared features, 1,124 stay
within the general `5.4e-13`, and 11 ill-conditioned ones depart by up to
`1.1156e-12`, `6.5974e-8` and `4.9245e-4` relative (`NONFINITE_FIT_GAP` in the
test).

**The tool.** `FeatureFinderCentroided` reads the mzML with the native reader,
which refuses non-finite binary values and scan start times; see
[TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT](TOPP_FEATURE_FINDER_CENTROIDED_SUPPORT.md),
native difference 13.

## Native differences

| Source behaviour | This port | Reason and evidence |
| --- | --- | --- |
| scores are float data arrays appended to each spectrum, replacing its existing float arrays | `ScoreArrays`: one flat `f32` array per score, outside the spectra | Only the algorithm reads them, and only debug mode writes them. The input spectra keep their arrays, and no per-spectrum allocation is needed. |
| `spectrumRanges().byMSLevel(1)` needs `updateRanges()` from the caller | ranges are computed on demand from the validated spectra, with `RangeBase`'s `std::min`/`std::max` semantics | No stale-range state exists, so the source's FAIMS "No ranges for this MS level" crash cannot occur. The source message "needs updated ranges" belongs to the peak-count check and is kept verbatim. |
| infinite and NaN retention times, m/z values, intensities and user-seed positions are read without a check | read as the Linux x86_64 Release build reads them; refused only at a NaN sort key whose order the port does not reproduce, at the endless profile merge of a NaN retention time and at the uncaught step-3.3.5 exception | See *Non-finite input*: 189 executed cases, 161 returned runs and 16 exceptions reproduced, 9 returned runs refused at a NaN sort key. The port refused every non-finite value before. |
| `mass_trace:min_spectra = 1` gives `min_spectra_ = 0`. Every trace score becomes 0/0 = NaN and every peak a local maximum; no overall score reaches a threshold, no seed is found and the run returns an empty map | the same: NaN trace scores, no seed, an empty map | The execution (B6 driver, `ffc1_min_spectra_1`) shows the source is *defined* here, so the port follows it (lead decision of 2026-09-15, `CPP-271`). B6 refused the configuration; that refusal is gone. Nothing later in the algorithm is reached, so the source's `size_t(-1)` delta buffer in `extendMassTrace_` stays unreachable; the port returns `Error::InvalidValue` if it ever is. |
| a zero or infinite intensity bin step makes `intensityScore_` convert `floor(NaN)` or `floor(inf)` to `UInt`, which is undefined | by default the Linux x86_64 Release build's outcome (every intensity score NaN, no seed, an empty map); `DegenerateBinStep::Refuse` refuses exactly the inputs whose seed loop reads those scores | Undefined behaviour of the `float`-to-`int` kind, whose Release outcome is measured, repeatable and explained by the emitted `cvttsd2si`; see *Degenerate intensity bins*. The port refused every zero-width range before, including short inputs, where the result does not depend on the scores. |
| `charge_low > charge_high + 1` wraps the `UInt` charge count and indexes past the score arrays | `Error::InvalidValue` from `Settings::charge_count` | Undefined behaviour. `charge_low == charge_high + 1` gives zero charges, as in the source. |
| `write_debug = true` writes `debug/` into the working directory, and `writeFeatureDebugInfo_` reads the undeclared `debug:pseudo_rt_shift`, which terminates the process at the first seed that reaches the fit | `debug::DebugOutput`, which the tool writes; `Error::Unsupported` at that seed under `PseudoRtShiftKey::Source` | Library code does not write files, and a safe port cannot end the process (see *Debug mode*). |
| a changed `abundance_12C` or `abundance_14N` builds the override from a default `IsotopeDistribution` that already holds `(0, 1)`; the patterns grow (FFC_1 window 0: 27 normalised bins) and FFC_1 with 12C = 90 % finds 0 seeds, 0 candidates and 0 features (C2 `ffap_ffc1_abundance_12C_90`) | the intended two-isotope distribution, which **does** find seeds and features | `CPP-247`, lead decision of 2026-09-15: follow the intent, not the defect, because the generator rejects the stray-peak construction and refusing a parameter the source accepts is worse. **This is the one place where the port's features differ from the executed C++ by design.** `AbundanceOverride::Refuse` is the opt-in for a caller that must not diverge. |
| the overall score is `std::pow(float, float)`, the platform `powf` | `libm::pow` in `f64`, rounded once to `f32` | The port's value is the correctly rounded power (60-digit decimal check) and the same on every machine. The reference build's glibc 2.39 `powf` (its FMA variant on the AMD EPYC 7763 capture host) is one binary32 step away on 8 of the 30,840 retained overall scores (2 of 3,084 on FFC_1); the macOS arm64 product SDK's Apple `powf` was on 99. `libm::powf` differed from Apple's on 226 of 3,084. No seed list changes. `tests/data/feature_finder_picked/overall_rounding.tsv` and the `rounding` rows of `degenerate_stage.tsv.gz` list every difference. |
| `std::sort` of seeds with equal `f32` intensity leaves an order the standard does not specify | the order the Linux Release build's libstdc++ introsort gives (`source_sort::source_sort_reversed_by`) | Deterministic, and equal to the reference platform (driver mode `sort`, 130 inputs). |
| `std::sort` of a bin's intensities: `-0.0` and `+0.0` compare equal | `total_cmp`: `-0.0` first | Only a quantile's zero sign can differ, and only for inputs with both signed zeros; the tool path filters non-positive intensities. |
| progress, `Found N seeds for charge c.` and `Found N feature candidates for charge c.` go to `std::cout`; the overlap count, the abort reasons, the feature count and the apex warning to `OPENMS_LOG_INFO`/`WARN` | `SeedStage::log` and `RunOutput::log`, in the source's order | Library code never prints. The candidate line is inserted directly after its charge's seed line, so the two `std::cout` lines are adjacent as in the source, even though this port computes every charge's seeds first. The bare newlines the source logs around the abort block (`FeatureFinderAlgorithmPicked.cpp:1019` and `1026`) are emitted as empty log entries, so a caller that prints the log line by line — `FeatureFinderCentroided` does — reproduces the executed C++ stdout block exactly. |
| steps 3.1 to 3.3 run per charge | steps 3.1 and 3.2 run for every charge first | Step 3.3 reads only its own charge's arrays, so the arrays and seeds are identical. |
| user seeds are a copied `FeatureMap` sorted with `std::sort` | positions only, sorted stably by `total_cmp`; a NaN m/z among two different other m/z values is `Error::InvalidValue` | Equal m/z values, and a NaN among copies of one value, give the same search result in any order; with two different other values the source sort is undefined. Infinite positions and a NaN retention time are ordinary (executed, 7 features). |
| unbounded work | `Limits`: spectra, peaks, charges, bins per dimension, windows, pattern values, score bytes, work units | Checked before the allocation or computation each bounds. The FFC_1 workload is several orders of magnitude below every default. |
| `Math::pearsonCorrelationCoefficient` returns an infinity when its denominator underflows to 0 from non-zero deviations | NaN, which counts as 0 | Inherited from `src/math/statistic_functions.rs`. Unreachable at isotope intensity scales. |
| `Exception::UnableToFit` from `fitter->fit` would be thrown inside the `omp parallel for`, where nothing catches it, and would end the process | unreachable, so not reproduced; an error from the port's fit (its point, byte or work ceilings, or the refused NaN profile merge) ends the run with that error | Both throws of `TraceFitter::optimize_` are unreachable from the seed loop: a fitted candidate has at least two traces of which at most one has fewer than three peaks, so at least 4 residuals for at most 4 parameters (`TraceFitter.cpp:111`), and Eigen returns `ImproperInputParameters` (`:129`) only for `maxfev <= 0`, which the `fit:max_iterations` restriction excludes. The residual count is `int` (`GaussTraceFitter.cpp:140`, `EGHTraceFitter.cpp:29`), so more than `INT_MAX` peaks would also throw at `:111`, but the solver's `MAX_POINTS` ceiling refuses such traces first. The argument is written out at `FittedModel::fit`. The port used to turn any fit error into an abort reason, which hid its own ceilings. |
| `extendMassTraces_` dereferences `pattern.spectrum[0]` when the pattern matched no peak | `Error::InvalidValue` | Undefined behaviour. Unreachable from `run_`, where the pattern always contains the seed; reachable through the public function, which the test exercises. |
| `traces[traces.max_trace]` is indexed before the fit without a range check | `Error::InvalidValue` | Undefined behaviour when `max_trace` is stale. The one branch that could make it stale (`traces.clear()` for a trace before `max_trace`) is unreachable, because `max_trace` is still 0 at every index that could satisfy `p < max_trace`. |
| `setWidth` stores any FWHM, including a NaN produced by a non-finite fit that the quality checks let through (every comparison against a NaN is false) | `Error::InvalidValue` from `BaseFeature::set_width` | The kernel setter validates. No executed configuration produces a non-finite fit; the 1,201 features of the 170 non-finite captures that returned, 140 of those captures with an infinite intensity, are all finite. |
| `f2.getCharge() % f1.getCharge()` divides by zero for a zero charge | `Error::InvalidValue` at that pair, only where the remainder the source evaluates would trap (also `INT_MIN % -1`) | The algorithm's own charges are at least 1. A caller's feature of charge 0 traps in the source (executed: SIGFPE, 2 of 2); same-charge pairs and non-overlapping charge-0 features are processed (executed). |
| a hull with no point has the default `DBoundingBox` `[DBL_MAX, -DBL_MAX]`, whose `width()` is negative infinity | the same arithmetic (`resolution.rs` `SourceBox`) | Only a caller's map can hold such a hull (executed: overlap cases). |
| `plot_nr` is assigned in an OpenMP critical section, so its value depends on the schedule | assigned in seed order | It is overwritten by the feature number for every feature that survives. The debug file names use it, and with one thread the source's value is this seed-order number. |
| `aborts_[reason]++` runs inside the parallel region without synchronisation | aggregated serially in seed order | A data race (candidate 5): with two or more threads aborting seeds the source has no reproducible result, so the behaviour at more than one thread cannot be tested against the C++ and is not. The port's counts are exact and independent of the thread count; every executed comparison of `aborts_` (C2, `degenerate_stage`) records the library's map at one thread only, and the multi-thread runs compare everything else. |
| the seed loop is an OpenMP `parallel for` with four named critical sections | `concept::parallel::map_collect` over the seed indices with `Options::threads` | `map_collect` returns results in input order, so no critical section is needed and the output is bit-identical at 1, 2 and 8 threads (`the_seed_loop_is_bit_identical_across_thread_counts`). The source's results are schedule-independent for the same reason, its `tmp_feature_map` being a `std::map` keyed by seed index; the C++ oracle gave identical output at 1, 2 and 8 threads. |
| the containment pass scans a growing `std::vector<Size>` of swallowed seeds | a `BTreeSet` | A pure membership test; duplicates in the source's vector change nothing. |
| `FeatureMap::sortByMZ` and `sortByIntensity` use `std::sort`, which is unstable | the Linux Release build's introsort order (`source_sort::source_sort_by`), NaN keys included; refused only at a step that would read outside the map | Equal to the reference platform (driver modes `sort`, `overlap`). |
| `setMetaValue("spectrum_index", Size)` stores an unsigned value | `i64`, refused above `i64::MAX` | The source's `DataValue` narrows to a signed integer anyway; the featureXML writer writes the same digits. |
| `CoarseIsotopePatternGenerator` iterates elements in heap-address order, which varies between runs | ascending atomic number | Inherited from B2. It is the majority order (198 of 200 runs), which every retained execution used; all windows match. |
| `MSExperiment::sortSpectra` and `sortChromatograms` with `std::sort` | module-local sorts that follow them (`std::is_sorted` first per spectrum, stable within), refusing only NaN keys (see *Non-finite input*) | The kernel's sorts refuse every non-finite value. Spectra with equal retention times and chromatograms with equal product m/z keep their input order, where `std::sort` leaves it unspecified; a spectrum with inconsistent attached records is an error. |
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

Every float of the **seed stage** is compared bit for bit on every platform;
the only exception is the documented overall-score rounding, whose correctly
rounded values the test asserts instead. The **feature stage** compares counts,
identities and orders exactly, and fitted quantities bit for bit on Linux
x86_64 with glibc except in the EGH configuration; the measured agreement is in
the next section.

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

**Linux x86_64 (glibc), measured on dax.** Since lane B3b the port's
Levenberg-Marquardt solver follows the Release build's Eigen kernels, and the
Gaussian fit calls the platform `exp` and `log`, as the source does. Every
compared value of the five Gaussian configurations is bit-identical to the
Release capture, including seeds 11 and 12 of
`classtest_9247_tight_pattern`, which departed by up to `2.3e-3` from the
product SDK and were the `KNOWN_FIT_GAP` of the earlier comparison. In the EGH
configuration (`ffc1_asymmetric`) 74 compared values depart, by at most
`2.3038e-12` relative (seed 24's lower retention-time bound; the fitted `tau`
of 16 seeds by up to `1.4e-13`; the features' `EGH_tau`, `EGH_sigma`,
`score_fit` and `score_correlation` by up to `7.9e-16`, while their retention
times, m/z values, intensities, qualities and widths are bit-identical):
`EGHTraceFitter` calls
the `libm` crate's `exp`, `log` and `atan` where the source calls glibc's
([EGH support](EGH_TRACE_FITTER_SUPPORT.md)). The test therefore compares
exactly on Linux x86_64 with glibc and bounds the EGH configuration by
`EGH_LIBM_GAP = 2.4e-12`. The exact comparison assumes a CPU with FMA, as both
measured hosts have: glibc selects FMA variants of `exp` and `log` there.

Counts, charges, labels, `num_of_datapoints`, hull counts, hull point counts,
hull point coordinates, subordinate counts, abort reasons and abort counts are
compared exactly and agree everywhere. Every isotope-fit score, isotope-pattern
intensity and m/z score and every mass trace (peak identities, theoretical
intensities, baseline) is bit-identical on every platform.

**macOS arm64, a platform note.** Against the same Linux capture, Apple's `exp`
in the Gaussian fit moves the last bits: every Gaussian value stays within
`5.355e-13` relative, except the two ill-conditioned fits of seeds 11 and 12 of
`classtest_9247_tight_pattern`, which depart by up to `1.07e-3` (area; height
`5.6e-4`, sigma and FWHM `5.1e-4`, upper bound `3.7e-5`, centre `2.2e-5`). Both
seeds are rejected by `checkFeatureQuality_` in the executed C++ and on macOS,
with the same reason, so no feature changes; the test asserts that a seed with
a platform gap never becomes a feature. The EGH configuration departs by the
same `2.3038e-12` as on Linux. On other platforms, which were not measured, the
test keeps the work package's `1e-9` contract and the `2.3e-3` bound measured
for the two seeds before the Linux capture existed.

**The intended abundance override (adapted).** The library cannot compute
the override the source intends, so `../oracle/ffap-sem-completion/drivers/intended_abundance.cpp`
recomputes step 2.5 with a cleared override distribution, assigns the windows
to the protected `isotope_distributions_` and replays steps 3.1 to 4 with the
library's protected functions (two repetitions at one and four threads,
identical). For FFC_1 with `abundance_12C` 90 and 99 and `abundance_14N` 95 the
port's default reproduces every window bit for bit, the seeds with their
pattern and overall scores (18, 25 and 13), the candidates, the abort reasons
and the features (1, 8 and 2), under the same platform bounds
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
  - the refusals: charge wrap, a NaN user-seed m/z among
    different ones, an unsorted spectrum with a NaN m/z, abundances under
    `AbundanceOverride::Refuse`, and degenerate bin steps under
    `DegenerateBinStep::Refuse` only when the seed loop reads them;
  - a single scan with a NaN intensity: sorted alone, scored with the zero
    steps' default NaN, and never reaching the seed loop (derived);
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
  --test isotopes_source_precision --test mass_trace --test mass_trace_detection \
  --test topp_feature_finder_centroided --test feature_finder_picked_instrumentation \
  --test progress_logger
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
   `double` and rounding once gives the correctly rounded value.
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
  - `charge_low > charge_high + 1` (the `UInt` charge count wraps and the score
    arrays are indexed past their end: an out-of-bounds access, refused where
    the count is computed).
  - Non-finite input is read as the Release build reads it (189 executed
    cases). What remains refused is listed in *Non-finite input*: a NaN
    sort key whose order the port does not reproduce at these sites
    (undefined, or an observable unspecified order; lifted once `source_sort`
    orders these sorts), the endless NaN
    profile merge, and the step-3.3.5 exception that terminates the source.
  - The `Limits` ceilings (bounded work).
  - `AbundanceOverride::Refuse` is an opt-out; the default computes the
    intended override, the one designed difference (`CPP-247`), now pinned
    against an adapted Release replay.
- **Unreachable source behaviour.** `Exception::UnableToFit` cannot be thrown
  from the seed loop (argument at `FittedModel::fit`, including the `int`
  residual count), and the `aborts_` data race at more than one thread has no
  reproducible C++ result, so neither is tested against the C++.
- Rust files: `src/analysis/feature_finder_picked/algorithm.rs`, `scoring.rs`,
  `seeds.rs`, `extension.rs`, `fitting.rs`, `resolution.rs`, `defs.rs`,
  `instance.rs`, `debug.rs` and `source_sort.rs`. Tests:
  `tests/feature_finder_picked_seeds.rs`, `tests/feature_finder_picked.rs`,
  `tests/feature_finder_picked_instrumentation.rs` and
  `tests/topp_feature_finder_centroided.rs`.
- The ledger scope's "steps 3.3 to 5" should read "steps 3.3 and 4": the
  source's last step is step 4 (`.cpp:860`).
- `docs/doc-coverage.json` needs `--write` for `defs.rs`, which is at 100 %;
  the floor was not re-recorded here because the file is the integrator's.
