# FeatureFinderCentroided: supported source subset

The TOPP tool that detects two-dimensional features in centroided LC-MS data.
Ported by package **C5-FFC-WRAPPER** of the early TOPP bundle, and reconciled
with the completed algorithm (package **B7-FFAP-FEATURES**) afterwards, from
`OpenMS4-topp/src/FeatureFinderCentroided.cpp` at topp `174b576`, with the
framework pieces it uses read at cli `c19e494` and the library pieces at core
`bc9cc12`. The Rust tool is `openms::cli::tools::FeatureFinderCentroided`
(`src/cli/tools/feature_finder_centroided.rs`); `src/bin/FeatureFinderCentroided.rs`
is the three-line executable, and `tests/topp_feature_finder_centroided.rs` runs
that executable. Hashes, oracle case names and the derivation rules of the
synthetic inputs are in
[`tests/data/topp_feature_finder_centroided_provenance.json`](../tests/data/topp_feature_finder_centroided_provenance.json).
The framework behaviour this tool inherits is documented in
[TOPP_CLI_SUPPORT](TOPP_CLI_SUPPORT.md), and the algorithm it wraps in
[FEATURE_FINDER_PICKED_SUPPORT](FEATURE_FINDER_PICKED_SUPPORT.md).

## What runs today

| Stage | State |
|---|---|
| Registration, `-write_ini`, INI merge and validation | complete; `-write_ini` is compared with the executed C++ file line by line with exact numbers and as a decoded parameter tree. `--help` is the framework's rendering of the same registration and has no oracle snapshot here |
| Loading `-in` (mzML, MS1, executed intensity range), `-seeds` (featureXML) | complete |
| Empty-input, per-peak ion-mobility and profile branches | complete, in source order |
| FAIMS input | **complete** (package B11): the split by compensation voltage, one algorithm run and seed filter per voltage, the `FAIMS_CV` annotation and the cross-voltage merge of `-faims_merge_features`, with the two defects of the source merge corrected (native difference 1) |
| The picked algorithm | complete (package B7); `-threads` reaches its seed loop |
| Primary MS run path, unique ids, `QUANTITATION` processing record, hull and subordinate clean-up, featureXML store | complete; reached by every accepted run and additionally tested on its own |
| `TOPP_FeatureFinderCentroided_1`, `-seeds`, `feature:rt_shape asymmetric`, `-debug 5`, `-threads 0/1/2/4/8` | run end to end and measured against the C++ Release build (see below); the `1e-9` acceptance against the C1 oracle outputs and the `-algorithm:fit:max_iterations` sweep are package B10's |

Every branch of `main_` therefore runs. The one reason the tool stayed
**partial** in the ledger — the FAIMS refusal — is gone, and the ledger scope
change is requested with this package. `TOPP_FeatureFinderCentroided_1`
produces the retained expectation.

## How close the output is to the executed C++

`TOPP_FeatureFinderCentroided_1` was run through this port and through
`/ceph/ibmi/abi/oliver/opt/openms4-release-bc9cc12-c19e494-174b576/bin/FeatureFinderCentroided`
on the same input and INI, and the two featureXML files compared decoded
(decision D6, generated identifiers excluded). Both give eight features of
charge 2 with 30 four-point hulls, 120 hull points, no subordinates, the same
`UserParam` key set, the same `spectra_data` and one `Quantitation` processing
record, and the same stdout (`25 seeds`, `8 feature candidates`, `Removed 0
overlapping features.`, `Invalid fit: Fitted model is bigger than 'max_rt_span':
1 times`, `8 features found.`).

| Field | Port vs C++ Release | C++ Debug vs C++ Release |
|---|---|---|
| convex-hull `rt` and `mz` (all 120 points) | bit-identical | bit-identical |
| feature `mz` | bit-identical | bit-identical |
| `intensity`, `FWHM` | identical as `f32` | identical as `f32` |
| `label`, `num_of_datapoints`, `spectrum_index`, `spectrum_native_id` | identical | identical |
| feature `rt` | `5.5e-13` relative | `2.2e-13` |
| `score_fit` | `2.2e-10` relative | `9.1e-11` |
| `score_correlation` | `7.7e-12` relative | `3.1e-12` |
| `overallquality` | agrees to the six decimals the C++ writer prints | — |

The residue is the last bits of the Levenberg-Marquardt fit, of the same order
as the C++ build's own Debug-to-Release spread on the same fields. These tool
measurements date from `4d53a7e`, before lane B3b made the fit follow the
Release build's Eigen kernels; the algorithm's features for FFC_1 are now
bit-identical to the Release build on Linux x86_64 with glibc
(`every_configuration_matches_the_executed_library` in
`tests/feature_finder_picked.rs`), and the tool was not re-measured here. The other
wrapper modes were measured the same way, on the same Linux node, against the
same Release build:

| Mode | Counts | Worst gap |
|---|---|---|
| `-seeds` with the retained output | `24 seeds`, `9 candidates`, `Removed 1 overlapping features.`, `Could not extend seed: 1 times`, 8 features, 30 hulls, 120 points | as FFC_1: `5.5e-13` rt, `2.2e-10` `score_fit`, `7.7e-12` `score_correlation` |
| `-algorithm:feature:rt_shape asymmetric` | 8 features, 30 hulls, 120 points, `EGH_height`/`EGH_sigma`/`EGH_tau` present | `rt` and `mz` bit-identical; `EGH_tau` `9.9e-16`, `EGH_sigma` `1.3e-16`, `score_fit` `6.2e-16`, `score_correlation` `5.6e-16` |
| `-debug 5` | 8 features, 30 hulls, **1054** points (no bounding-box reduction) | as FFC_1 |

`-threads 0`, `1`, `2`, `4` and `8` write byte-identical output, unique ids
included. The port's own two builds agree to one ulp: the macOS arm64 and Linux
x86_64 `FeatureFinderCentroided` differ by at most `1.1e-16` relative on the
FFC_1 output, so the numbers above are not a platform artefact.

## API mapping

Every member of the source class, and every framework call it makes.

| C++ | Rust |
|---|---|
| `class TOPPFeatureFinderCentroided : public TOPPBase` | `pub struct FeatureFinderCentroided` with `impl Tool` |
| `TOPPBase("FeatureFinderCentroided", "Detects two-dimensional features in LC-MS data.", true, {Citation…})` | `Tool::NAME`, `Tool::DESCRIPTION`; the `official` flag and the two citations are not ported (the framework ports neither) |
| `registerOptionsAndFlags_()` | `Tool::register` |
| `registerInputFile_("in", "<file>", "", "input file")` + `setValidFormats_("in", {"mzML"[, "raw"]})` | `ToolSpec::register_input_file` + `set_valid_formats`; the `raw` entry exists only in a `WITH_THERMO_RAW` build and is not ported |
| `registerOutputFile_("out", …)` + `setValidFormats_("out", {"featureXML"})` | `ToolSpec::register_output_file` + `set_valid_formats` |
| `registerInputFile_("seeds", …, false)` + `setValidFormats_("seeds", {"featureXML"})` | as above, `required = false` |
| `registerStringOption_("faims_merge_features", "<true/false>", "true", …, false)` + `setValidStrings_` | `ToolSpec::register_string_option` + `set_valid_strings`; read before the merge |
| `addEmptyLine_()` (twice) | `ToolSpec::add_empty_line` |
| `registerSubsection_("algorithm", "Algorithm section")` | `ToolSpec::register_subsection` |
| `getSubsectionDefaults_(const std::string&)` | `Tool::subsection_defaults`, returning `crate::analysis::feature_finder_picked::algorithm::default_parameters` for any name, as the source ignores the name |
| `main_(int, const char**)` | `Tool::run_io`; `Tool::run` forwards the process streams |
| `getStringOption_("in" / "out" / "seeds")` | `ToolContext::string` |
| `PeakFileOptions` with `setMSLevels({1})` and `setIntensityRange({min(), maxPositive()})` | `FeatureFinderCentroided::peak_file_options` |
| `FileHandler::loadExperiment(in, exp, {MZML, RAW}, log_type_)` | `FileHandler::load_experiment_with_read_options` with those options, the same allowed types and the source-compatible mzML reader options (decision D10); the progress log type is the framework's and no loader reports progress |
| `exp.updateRanges()` | nothing: native ranges are computed on demand |
| `throw FileEmpty("Error: No MS1 spectra in input file.")` | `FeatureFinderCentroided::NO_MS1_SPECTRA_MESSAGE` and `ExitCode::InputFileEmpty` |
| `IMTypes::determineIMFormat(spec) == IMFormat::IM_PEAK` | `ImTypes::determine_im_format` over every spectrum |
| `imPeakTypeToString(spec.getIMPeakType())` | `im_peak_type_to_string(IonMobilityPeakType::Profile)` in `FeatureFinderCentroided::im_peak_message` (see the native differences) |
| `return INCOMPATIBLE_INPUT_DATA` | `ExitCode::IncompatibleInputData` |
| `exp[0].getType()` and `getFlag_("force")` | `experiment.spectra[0].spectrum_type` and `ToolContext::force` |
| `throw IllegalArgument("Error: Profile data provided …")` | `FeatureFinderCentroided::PROFILE_DATA_MESSAGE` and `ExitCode::UnknownError`, the code `TOPPBase` gives that exception |
| `FileHandler().loadFeatures(seeds, map, {FEATUREXML})` | `FileHandler::load_feature_map` |
| `getParam_().copy("algorithm:", true)` | `ToolContext::subsection("algorithm")` |
| `writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3)` | not ported: the framework ports no debug log or `-log` file |
| `IMDataConverter::splitByFAIMSCV(std::move(exp))` | `ImDataConverter::split_by_faims_cv` (B8), whose `FaimsSplit::messages` are written to the two log-stream caches by their level. The source's `std::pair<double, MSExperiment>` key is `FaimsGroupKey`, whose `NotFaims` variant names the source's NaN |
| `const bool has_faims = faims_groups.size() > 1 \|\| !std::isnan(faims_groups[0].first)` | `FaimsSplit::has_faims`, the same predicate on the named key |
| the per-group loop `for (auto& [group_cv, faims_group] : faims_groups)` | the loop over `FaimsSplit::groups`, ascending by voltage as the source's `std::map` is |
| `OPENMS_LOG_INFO << "Processing FAIMS CV group: " << group_cv << " V (" << faims_group.size() << " spectra)"` | `FeatureFinderCentroided::processing_group_message`, the voltage at the stream's default precision 6 |
| the seed filter (`258-281`), `Constants::UserParam::FAIMS_CV`, tolerance `0.01` | `FeatureFinderCentroided::seeds_of_group`, with `FaimsHelper::DEFAULT_CV_TOLERANCE` |
| `feat.setMetaValue(Constants::UserParam::FAIMS_CV, group_cv)` | the same meta value on each of the group's features, before they join the result |
| `OPENMS_LOG_INFO << "Combined " << features.size() << " features from all FAIMS CV groups."` | `FeatureFinderCentroided::combined_features_message` |
| `getStringOption_("faims_merge_features") == "true"`, `FeatureOverlapFilter::mergeFAIMSFeatures(features, 5.0, 0.05)` and its `"FAIMS feature merge: …"` line | `FeatureOverlapFilter::merge_faims_features_with_fidelity` (B9) with `FAIMS_MERGE_MAX_RT_DIFF`, `FAIMS_MERGE_MAX_MZ_DIFF` and `FeatureFinderCentroided::FAIMS_MERGE_FIDELITY`, then `FeatureFinderCentroided::faims_merge_message`. The fidelity is the tool's one designed difference here (native difference 1) |
| `OPENMS_LOG_INFO << "Not FAIMS compensation voltages found …"` (in `splitByFAIMSCV`) | `FeatureFinderCentroided::NO_FAIMS_MESSAGE` on the output stream |
| `OPENMS_LOG_INFO << "FAIMS data detected with N compensation voltage(s)."` | `FeatureFinderCentroided::faims_detected_message` |
| `FeatureFinderAlgorithmPicked ff; ff.run(std::move(group), features_cv, feafi_param, seeds_cv)` | `crate::analysis::feature_finder_picked::instance::FeatureFinderAlgorithmPicked::with_options(…).run(…)` into an empty map, then `take_debug_output` and `report`, in `run_group`; one fresh object per group, as the loop body creates one |
| the files `ff.run` writes with `-algorithm:write_debug` (`debug/log.txt`, `debug/seeds_<charge>.featureXML`, `debug/abort_reasons.featureXML`, `debug/input.mzML`, `debug/features/*`) and the `FeatureXMLHandler::store()` lines they log | `write_debug_log`, `store_debug_features` and `FeatureFinderCentroided::run_io`, from the algorithm's `DebugOutput`, in the source's order (see *Debug mode*) |
| `features.setPrimaryMSRunPath({"file://" + File::basename(in)})` under `-test`, else `{in}` | `FeatureFinderCentroided::finish_features` step 1 (`FeatureMap::set_primary_ms_run_path`, `file::basename`) |
| `features.ensureUniqueId(); features.applyMemberFunction(&UniqueIdInterface::setUniqueId)` | step 2 (`HasUniqueId::ensure_unique_id` on the map's id, then `FeatureMap::for_each_unique_id` with one `ToolContext::unique_id_generator`) |
| the `debug_level_ > 10` metadata listing | step 3 |
| `addDataProcessing_(features, getProcessingInfo_(DataProcessing::QUANTITATION))` | step 4 (`ToolContext::processing_info` + `ToolContext::add_data_processing`) |
| `debug_level_ < 5`: `ft.getConvexHull().expandToBoundingBox()`, every `getConvexHulls()[i].expandToBoundingBox()`, `ft.getSubordinates().clear()` | step 5, without the overall hull (see the native differences) |
| `map_file.storeFeatures(out, features, {FEATUREXML})` | `FileHandler::store_feature_map` |
| `return EXECUTION_OK` | `ExitCode::ExecutionOk` |
| `int main(int argc, const char** argv)` | `src/bin/FeatureFinderCentroided.rs` |

## Exit codes and diagnostics

Every row was executed against the product SDK; the case name is the oracle
case (C1 = `../oracle/topp-early-bundle`, C5 = `../oracle/ffc-wrapper-c5`). The
two huge-m/z rows were executed against the Release build
(F2 = `../oracle/ffap-complete-fix2`, two runs each, identical apart from the
timing text).

| Input | Exit | Diagnostic | Oracle case |
|---|---|---|---|
| Registered workflow `TOPP_FeatureFinderCentroided_1` | 0 | — (eight features, see above) | C1 `TOPP_FeatureFinderCentroided_1` |
| `-write_ini` (with and without `-test`) | 0 | — | C1 `FFC_write_ini`, `FFC_write_ini_test` |
| MS2 spectra only | 4 | `Error: File empty (the file 'Error: No MS1 spectra in input file.' is empty)` | C1 `FFC_ms2_only` |
| Per-peak ion mobility, with or without a unit | 11 | `Error: Input contains per-peak ion mobility data (IM_PEAK, im_profile) …` | C1 `FFC_im_peak_with_units`, `FFC_im_peak_without_units` |
| Per-peak ion mobility **and** profile data | 11, the ion-mobility message | the check order of `main_` | C5 `c5_im_peak_profile_noforce` |
| Ion-mobility arrays on MS2 spectra only | 0, empty feature map | — | C1 `FFC_im_arrays_ms2_only` (its Debug exit 8 is a precondition, `debug_only`); the C++ Release build exits 0 with `0 features found.`, as this port does |
| First spectrum stored profile, no `-force` | 8 | `Error: Unexpected internal error (Error: Profile data provided but centroided spectra expected. …)` | C1 `FFC_profile_noforce`, `FFC_FileFilter_44_noforce`, C5 `c5_first_profile_only` |
| Two MS1 spectra (`FileFilter_44_input.mzML` with `-force`) | 0, empty feature map, no seed for charges 1 to 4 | — | C1 `FFC_FileFilter_44_force` (its Debug exit 8 is a precondition, `debug_only`); the C++ Release build exits 0 with the same lines and map, as this port does; see native difference 2 |
| First spectrum profile, `-force` | 0, the FFC_1 features | — | C1 `FFC_profile_force` |
| Profile term before `MS:1000525`, no `-force` | 0, the FFC_1 features: the reader resets the type | — | C1 `FFC_profile_then_spectrum_representation` |
| Every spectrum but the first stored profile | 0, the FFC_1 features: only `exp[0]` is checked | — | C5 `c5_later_profile_only` |
| Every MS1 peak below the intensity range | 8 | `Error: Unexpected internal error (FeatureFinder needs updated ranges on input map. Aborting.)` | C5 `c5_negative_intensities` |
| `-seeds` that is not featureXML | 6 | `Input file '…' has invalid format 'mzML'. Valid formats are: 'featureXML'.` | C5 `c5_seeds_not_featurexml` |
| FAIMS input with an unreadable `-seeds` file | 3 | `Error: Unable to read file (…)`, before any FAIMS message | C5 `c5_faims_corrupt_seeds` |
| Any FAIMS input (one, two or three voltages, a voltage on half the spectra, `-faims_merge_features false`, the two upstream FAIMS fixtures) | C++ 8; this port 0 | C++ `Error: Unexpected internal error (the value '1' was used but is not valid; No ranges for this MS level)` after the first `Processing FAIMS CV group:` line and no output (`CPP-278`); this port writes the features of every voltage group, merged unless `-faims_merge_features false` | C1 `FFC_faims_*`, C5 `c5_faims_partial_cv`, B11 `faims_*` (`../oracle/b11-faims`, seven inputs, all rc 8); native difference 1 |
| A FAIMS input whose skipped spectra carry no voltage | 0 | `Skipping spectrum without FAIMS CV (no prior FAIMS CV context or unexpected layout).` on **stderr**, folded by the warn stream's cache into one line and `<…> occurred 56 times`, as the executed C++ wrote it | C5 `c5_faims_partial_cv`, B11 `faims_partial_cv` |
| FAIMS **profile** input without `-force` | 8, the profile message | the profile check precedes the split | C1 `FFC_faims_interleaved_noforce` |
| `-algorithm:feature:rt_shape bogus` | 6 | `Invalid string parameter value 'bogus' … Valid values are: 'symmetric,asymmetric'.` | C1 `FFC_invalid_rt_shape` |
| `-out` without an extension | 0, the FFC_1 features written into it | — | C1 `FFC_out_no_extension` |
| FFC_1 with its last m/z `1e19`, `-algorithm:write_debug` (step 2.5 needs `ceil(1e19 * 2 / 100) + 1 = 2e17 + 1` isotope windows at the INI's charge 2 and width 100, more than `vector::max_size()`) | 12 | `Unable to initialize or run FeatureFinderCentroided: vector::_M_default_append`: the `std::length_error` reaches TOPPBase's outer `std::exception` handler (`TOPPBase.cpp:519-522`); `debug/features` and a 40-byte `debug/log.txt` are left | F2 `tool_1e19`; `a_debug_run_beyond_the_isotope_window_limit_exits_as_the_release_build` |
| the same with its last m/z `2e18` (`4e16 + 1` windows) | C++ 12; this port 8 | C++ `Unable to initialize or run FeatureFinderCentroided: std::bad_alloc`; this port `Error: Unexpected internal error (… exceed the limit of 1000000)`, the same debug files | F2 `tool_2e18`; native difference 14 |
| FFC_1 INI with `feature:reported_mz average` and the trace, seed, feature and isotope-fit thresholds at 0 (a seed whose best isotope pattern stayed empty reaches `extendMassTraces_`) | C++ SIGSEGV, shell status 139; this port 8 | C++ nothing on stderr and no output, the console lines up to the FAIMS line; this port `Error: Unexpected internal error (FeatureFinderAlgorithmPicked seed extension: the isotope pattern matched no peak; the source reads its first entry here)`, no output | `avg0` (`../oracle/ffap-complete-fix3`); `a_run_that_reaches_an_empty_best_pattern_is_refused_where_the_release_build_crashes`; native difference 15 |
| FFC_1 with every retention time scaled by `1e36` or `1e39` in the text (features of infinite intensity, and of infinite width at `1e39`) | C++ 0; this port 3 | C++ none, a featureXML with `inf` values; this port `Error: Unable to read file (parse error on line 0: nonfinite feature value)`, no output | `rt_e36`, `rt_e39`, `rt_e39_egh` (`../oracle/ffap-complete-fix4`); `infinite_feature_values_are_refused_by_the_featurexml_writer`; native difference 16 |

No branch writes `-out` before the store, and every refusal above was checked
to leave no output file.

## Debug mode

`-algorithm:write_debug` is a `true`/`false` string parameter, which TOPPBase
registers as a flag: it takes no value, and `-algorithm:write_debug true` ends
in exit 6 with `Trailing arguments after flag` (case x0, matched). With the
flag the algorithm returns its debug output (see
[FEATURE_FINDER_PICKED_SUPPORT](FEATURE_FINDER_PICKED_SUPPORT.md), *Debug
mode*), and the tool writes it where the source writes it: it creates
`debug/features` in the working directory, writes `debug/log.txt`, and stores
the seed maps, the abort map and the input with the scores at the points of the
run where the source stores them, with the `FeatureXMLHandler::store():  found N
invalid unique ids` lines they log. The console output passes through a model
of OpenMS's log-stream line cache (`LogStream.cpp`), so repeated lines are
folded into `<line> occurred N times` as in the executed output. The abort
map's unique id is drawn from the `-test` generator before the output's ids,
as in the source, so the output ids match too.

Executed against the Release `FeatureFinderCentroided`
(`../oracle/ffap-instr-completion`, `tool_cases.py`, `OMP_NUM_THREADS=1`,
three repetitions unless stated):

| Case | C++ Release | This port | Compared |
| --- | --- | --- | --- |
| a1: FFC_1 INI, `mass_trace:min_spectra 1` (no seed) | exit 0 | exit 0 | the console lines from the FAIMS line on (without the `took` line), `-out` with its unique id, and every debug file: `log.txt` byte for byte, featureXML and mzML decoded (D6), the abort map's unique id |
| a2: FFC_1 INI, `feature:min_isotope_fit 1.0` (every seed aborts) | exit 0 | exit 0 | as a1 |
| a3: `FileConverter_31_output.mzML`, `-force`, default parameters | exit 0 | exit 0 | as a1, four seed maps and the `occurred 4 times` line |
| a4: FFC_1 input, default parameters with `feature:min_isotope_fit 1.0` (reaches the fit); 2 repetitions | SIGABRT, shell status 134 | the b1 path | at the library level (`feature_finder_picked_instrumentation`): `log.txt` up to termination, `seeds_1.featureXML` |
| b1: FFC_1 INI (reaches the fit) | SIGABRT, shell status 134, OpenMS's fatal-exception block on stdout | exit 8, `Error: Unexpected internal error (the element 'debug:pseudo_rt_shift' could not be found)` | the console lines before the block, `log.txt` up to the last byte the file buffer had written, `seeds_2.featureXML`, `debug/features/` present, no abort map, no input, no `-out` |
| c1, c2: as b1 and a2 with `-threads 4` | as b1, a2 | the single-thread files | recorded, not compared: the executed logs differ between repetitions (data race in the source) |
| c3: as a1 with `-threads 4` (no seed, so nothing is written inside the parallel region; 2 repetitions) | exit 0, every file identical to a1's | exit 0, the single-thread files | as a1: `log.txt` byte for byte against c3's, the other files against a1's, which the executed c3 files equal byte for byte |
| `tool_1e19`, `tool_2e18` (`../oracle/ffap-complete-fix2`, 2 repetitions): FFC_1 INI and `feature:min_isotope_fit 1.0` on FFC_1 with its last m/z `1e19` or `2e18`, which step 2.5 cannot allocate | exit 12 (`length_error`, `bad_alloc`) | exit 12 for `1e19`; exit 8 for `2e18` (native difference 14) | the stdout block, stderr, `debug/features` present and empty, `debug/log.txt` (40 bytes, the first log line), no other debug file, no `-out` |

**Termination.** Every write_debug run in which a seed reaches the fit
terminates the C++ process, because the algorithm reads an undeclared
parameter inside its OpenMP region. A safe port cannot abort the process; the
tool writes everything the executed run wrote before it died and exits 8 with
the message TOPPBase gives that exception where it can catch it. The tool
always uses the source's key (`PseudoRtShiftKey::Source`); a library caller can
choose the declared key instead.

## The FAIMS closure

### What the merge is meant to do

The cross-voltage merge has no C++ oracle: the executed C++ tool never reaches
it, so no C++ build produces a merged FAIMS feature. What it is *meant* to
produce is nevertheless written down, in the parameter documentation of
`FeatureOverlapFilter::mergeFAIMSFeatures` and `mergeOverlappingFeatures`
(`FeatureOverlapFilter.h`), and the derivation is this, sentence by sentence:

1. *Merge FAIMS features that represent the same analyte detected at different
   CV values.* A cluster is **one analyte**, so it collapses to **one**
   feature. The source's result for three voltages is two features.
2. *Features are considered the same analyte if they have DIFFERENT FAIMS_CV
   values, are within `max_rt_diff` seconds in RT, within `max_mz_diff` Da in
   m/z, and have the same charge state.* Membership is a property of the pair
   of features, not of how many merges have happened; after a merge the
   survivor stands for the set of voltages in `merged_centroid_IMs`, so the
   test *different CV* becomes *a voltage the survivor does not yet stand for*.
   The source instead tests the survivor's `FAIMS_CV`, which its own callback
   has just removed, so it refuses every merge after the first.
3. *The feature with highest intensity is kept, and intensities are summed.*
   The kept feature is the cluster's maximum, which the descending-intensity
   sort puts first, and the sum is the analyte's total — so each member is
   counted **once**. The source offers a feature it has already removed to the
   next survivor, which adds it a second time.
4. *`merged_centroid_rts` / `merged_centroid_mzs` / `merged_centroid_IMs`:
   positions of all features that were merged; `FAIMS_merge_count`: number of
   FAIMS CV values that were merged.* The lists therefore grow to the size of
   the cluster, starting with the survivor's own values, and the count is the
   length of the voltage list.
5. *Features without `FAIMS_CV` are left unchanged*, and *non-FAIMS data: no
   merging occurs; single-CV FAIMS data: no merging occurs.* Unchanged here.

Two consequences that the wording does not state and that the port fixes by
the same reasoning: a survivor is never absorbed into another cluster (it
carries no `FAIMS_CV` any more, and it is the member of higher intensity,
which point 3 keeps), and the arithmetic is the source's — each merge stores
`double + double` into a `float`, so the running sum is narrowed to `f32` after
every step and the order in which a survivor absorbs its cluster can change the
last bit.

`FaimsMergeFidelity::Corrected` is exactly that merge; `FaimsMergeFidelity::Source`
is the executed one, kept for a caller who must not diverge, and tested against
the executed `c2_*` cases of `tests/feature_overlap_filter.rs`.

### The oracle, built from the parts

Because no C++ run produces the whole-tool answer, each compensation-voltage
group is written as its own single-voltage mzML — the spectra the split assigns
to that voltage, in input order, with the FAIMS cvParam removed so the C++ tool
takes its non-FAIMS path — and the C++ Release build is run on that file. The
port's features for that group must equal that run's under decision D6. That
pins the split, the per-group algorithm run, the group order and the annotation
against executed C++; only the merge is left, and it is pinned against the
derivation above with hand-derived numbers.

Executed on ibminode06 against
`openms4-release-bc9cc12-c19e494-174b576`, `OMP_NUM_THREADS=1`, `-test`
(`../oracle/b11-faims`: `derive.py`, `split_groups.py`, `run.sh`, `cases.sh`,
`cases2.sh`):

| Case | Input | C++ Release | This port |
|---|---|---|---|
| `faims_one_cv` | FFC_1 input with `-45` on every scan, FFC_1 INI | exit 8 after `FAIMS data detected with 1 compensation voltage(s).` and `Processing FAIMS CV group: -45 V (112 spectra)` | exit 0; the group is the whole input, so the features are `TOPP_FeatureFinderCentroided_1`'s, each with `FAIMS_CV` `-45`; the merge merges nothing (one voltage) |
| `faims_two_cv`, `faims_two_cv_nomerge` | `-45` and `-60` in turn, FFC_1 INI | exit 8 after the `-60` group line | exit 0; `-60` gives 3 features and `-45` gives 2, each equal to `group_m60` and `group_m45` below; merged, 3 features |
| `faims_three_cv` | `-45`, `-60`, `-70` in turn, FFC_1 INI with `mass_trace:min_spectra 5` | exit 8 after the `-70` group line | exit 0; 8 + 8 + 7 features equal to `group3s5_*`; merged, 10 features |
| `faims_partial_cv` | `-45` on every second scan, FFC_1 INI | exit 8, with the skip warning 57 times on stderr | exit 0, the same stderr; the group equals `group_m45` |
| `faims_test_data`, `faims_interleaved` | the two upstream fixtures (`CPP-240`: they spell the volt unit `UO:000218`) | exit 8 | exit 0, an empty map, as the group runs below |
| `group_m45`, `group_m60` | the two voltage groups of `faims_two_cv` as their own mzML, FFC_1 INI | exit 0, 2 and 3 features | the same, feature by feature (D6) |
| `group3s5_m45`, `group3s5_m60`, `group3s5_m70` | the three groups of `faims_three_cv`, `min_spectra 5` | exit 0, 8, 8 and 7 features | the same |
| `testdata_cvm65`, `interleaved_cvm45`, `interleaved_cvm60` | the upstream fixtures' groups | exit 0, `0 features found.` | the same |

The five group outputs are retained as fixtures
(`tests/data/topp_feature_finder_centroided/faims_group*.featureXML`); the
inputs are re-derived in the tests from `FeatureFinderCentroided_1_input.mzML`
by the recorded rule and pinned to the executed file by SHA-1.

### The merged numbers, hand-derived

Two voltages, from `group_m45` and `group_m60`; every pair is within 5 s and
0.05 Da and has charge 2, and the sum is the `f32` of the two `f64` intensities:

| analyte | −45 V | −60 V | survivor | merged |
|---|---|---|---|---|
| 4389.11 s, 648.257 Da | 45181.723 | 44601.215 | −45 | **89782.9375** |
| 4300.95 s, 651.760 Da | 35109.383 | 34660.066 | −45 | **69769.453125** |
| 4278.07 s, 653.770 Da | — | 19216.809 | −60 | **19216.80859375**, `FAIMS_CV` kept |

`5 -> 3 features (merged 2)`. Each survivor carries `merged_centroid_IMs`
`[-45, -60]`, `FAIMS_merge_count` 2, and the two retention times and m/z values
of its pair, and has lost its `FAIMS_CV`.

Three voltages with `min_spectra 5` are the case `CPP-283` answers: six of the
ten clusters hold **three** features. The source's merge would leave those six
as twelve features of 1900/1700 shape; the corrected merge gives one each,
`23 -> 10 features (merged 13)`. The ten survivors, their merged voltage lists
and their intensities are written out in
`three_faims_voltages_collapse_a_cluster_of_three_into_one_feature`. One
cluster there — 4278.16 s, 653.776 Da — is the only place in this package where
the `f32` running sum depends on the order in which the survivor absorbs its
two partners (58645.19921875 or 58645.203125). That order is the quadtree's
traversal order: deterministic for a given input, but not derivable by hand, so
the test accepts those two values and nothing else and says why.

## Reusing an instance

The source tool creates a fresh `FeatureFinderAlgorithmPicked` for each FAIMS
group and runs it into an empty map, so no state carries over between runs,
and so does this port. The object's behaviour across runs and with a caller's
non-empty map is ported and tested at the library level
([FEATURE_FINDER_PICKED_SUPPORT](FEATURE_FINDER_PICKED_SUPPORT.md), *Reusing an
instance*).

## Preserved source conventions

- **Check order.** Load, empty input, per-peak ion mobility, profile, seeds,
  `algorithm:` parameters, FAIMS, algorithm, annotation, clean-up, store. Three
  of these orders are observable and each is pinned by an executed case:
  ion mobility before profile, profile before FAIMS, seeds before FAIMS. The
  split is where the source performs it, so a FAIMS profile file without
  `-force` still exits 8 with the profile message and an unreadable `-seeds`
  file still exits 3 before any FAIMS line.
- **Only `exp[0]`.** The profile check reads the stored type of the first
  spectrum only, and never estimates it from the data. An mzML `MS:1000525
  spectrum representation` term after a profile term resets the stored type to
  unknown (`MzMLHandler.cpp:1634-1645`), so such a file is not refused.
- **The intensity filter keeps zero.** `main_` passes
  `std::numeric_limits<DRange<1>::PositionType>::min()`, and `DPosition<1>` has
  no `numeric_limits` specialisation, so the executed lower bound is `0` and not
  `DBL_MIN` as the source comment says. `peak_file_options` passes `0.0`; the
  half-open range still drops negative intensities. A2/A3 recorded the same
  reading (`MZML_MOBILITY_SUPPORT.md`), and `c5_negative_intensities` executes
  the consequence.
- **`FileEmpty` wording.** `TOPPBase` prints a `FileEmpty` exception as
  `Error: File empty (the file '<message>' is empty)`, so the source's message
  appears where a file name would; reproduced verbatim.
- **`IllegalArgument` is exit 8.** Both the profile refusal and the algorithm's
  own argument checks are `IllegalArgument`, which reaches `TOPPBase`'s
  `BaseException` arm, not its `InvalidParameter` arm. The wrapper maps the
  algorithm's `Error::InvalidValue` to `Error: Unexpected internal error
  (<message>)` with exit 8 for that reason.
- **`-test` outputs.** The primary MS run path keeps only the base name behind
  `file://`, the unique-id generator is seeded, and the processing record
  carries `version_string`, `1999-12-31 23:59:59` and `parameter: mode` =
  `test_mode`.
- **Unique ids.** `ensureUniqueId` draws one id for a map that has none, which
  the following `applyMemberFunction(setUniqueId)` immediately overwrites; the
  port makes the same two draws from one generator, so a later comparison with
  C++ ids is meaningless either way (the upstream comparison whitelists `id=`).
- **Clean-up below debug level 5.** Every mass-trace hull becomes its four-point
  bounding box and the subordinate features are dropped, *to reduce file size of
  feature files*; from level 5 on both are kept.

## Native differences

1. **FAIMS input is processed, with the two defects of the source merge
   corrected.** Decision D5 deferred the closure and the tool refused FAIMS
   input with exit 11; package B11 ships it. The corrected points, each named
   at the item it answers:

   - **`CPP-278`, the missing ranges, cannot arise here.** The source builds
     each voltage group with `addSpectrum` and never calls `updateRanges`, so
     `FeatureFinderAlgorithmPicked` throws `the value '1' was used but is not
     valid; No ranges for this MS level` on the first group and every FAIMS
     input exits 8. Re-executed for this package: seven FAIMS inputs through
     the Release build, all rc 8, each after printing `FAIMS data detected with
     N compensation voltage(s).` and the first `Processing FAIMS CV group:`
     line (`../oracle/b11-faims`). The native containers compute ranges on
     demand (`MSExperiment::spectrum_range_manager`), so a group has its own
     ranges the moment it holds spectra. There is no state to forget, nothing
     to emulate and no option: the defect is not reachable in this port. That
     is a property of the container port, not a choice made here.
   - **`CPP-282`, the merge that erases every feature.** `mergeFAIMSFeatures`
     records removal in an `unordered_set` keyed by `getUniqueId()`, and
     `FeatureFinderAlgorithmPicked` returns every feature with id 0, so the
     first merge marks id 0 removed, every later feature is skipped as a
     querier, and the closing `erase` drops all of them. Executed: oracle case
     `c2_uid0_wipe` of `../oracle/feature-overlap-filter`, three features in
     and none out. The tool draws a unique id for each feature from its
     generator **before** the merge, so removal keys on real ids. Those ids are
     overwritten immediately afterwards by the source's own
     `applyMemberFunction(&UniqueIdInterface::setUniqueId)`, so the only
     observable effect is that a FAIMS run consumes that many draws earlier
     than a non-FAIMS run; the ids of a featureXML are excluded from every
     comparison anyway (upstream whitelists `id=`). The corrected merge refuses
     a map whose FAIMS features share an id rather than silently erasing them.
   - **`CPP-283`, the double count and the survivor that stops absorbing.**
     `FeatureFinderCentroided::FAIMS_MERGE_FIDELITY` is
     `FaimsMergeFidelity::Corrected`, the library option B11 added beside the
     source-following `FaimsMergeFidelity::Source` (the pattern of
     `AbundanceOverride`, native difference 3). It skips a candidate already
     marked removed and lets a survivor keep absorbing voltages it does not yet
     stand for. See *What the merge is meant to do* below.

   Two further defects of the split belong to the ported library and are
   unchanged by this package: a NaN compensation voltage is refused rather than
   allowed to break an ordered set (`CPP-280`, native difference 4), and the
   chromatograms the source destroys are returned by `FaimsSplit` instead of
   dropped (`CPP-279`) — this tool loads MS level 1 only and uses none, so it
   discards them as the source does, and no output of this tool can show the
   difference. The information line of a non-FAIMS input keeps the source's
   wording, `Not FAIMS compensation voltages …` (`CPP-281`), because it is the
   line every executed run of this tool prints and every console comparison in
   this document rests on it.
2. **Degenerate intensity bins and short inputs follow the C++ Release
   build.** Where every MS1 spectrum has one retention time or every MS1 peak
   one m/z, the algorithm's intensity bin step is zero and the source converts
   `floor(NaN)` to `UInt` for every peak, which is undefined behaviour. The
   tool uses the library default `DegenerateBinStep::Source`, which reproduces
   what the C++ Release build computes there (every intensity score NaN, no
   seed), so the tool exits 0 with an empty map and the source's lines, as the
   executed Release tool does
   ([FEATURE_FINDER_PICKED_SUPPORT](FEATURE_FINDER_PICKED_SUPPORT.md),
   *Degenerate intensity bins*). Executed against
   `openms4-release-bc9cc12-c19e494-174b576` (`../oracle/ffap-sem-completion`,
   `tool_cases.py`, three repetitions each, identical), and asserted:

   | Case | Input | C++ Release and this port |
   |---|---|---|
   | `zero_rt`, `zero_rt_threads4`, `zero_rt_min_score_0` | FFC_1 with every `scan start time` 4114.53, FFC_1 INI; `-threads 4`; `seed:min_score` 0 | exit 0, `Found 0 seeds for charge 2.`, `Found 0 feature candidates …`, `0 features found.`, empty map |
   | `zero_mz`, `zero_mz_threads4`, `zero_mz_min_score_0` | FFC_1 with every m/z 500 | the same |
   | `zero_mz_control_min_score_0` | as `zero_mz` with one m/z 499, `seed:min_score` 0 | exit 0, `Found 735 seeds`, 0 candidates, `Could not find good enough isotope pattern containing the seed: 735 times`, empty map |
   | `filefilter_44_force` | `FileFilter_44_input.mzML`, `-force` | exit 0, no seed and no candidate for charges 1 to 4, empty map |
   | `fileconverter_31` | `FileConverter_31_output.mzML` | the same |

   `FileFilter_44_input.mzML` holds **two** MS1 spectra (at 0.273 s), and its
   result is fixed by its length, not by its zero extent: with the default
   `mass_trace:min_spectra` of 10 the seed loop is empty for any input of at
   most 10 scans. The C++ Debug build stops both short inputs with an
   `OPENMS_PRECONDITION` in `ProgressLogger::startProgress`
   (`ProgressLogger.cpp:235`), reached from `startProgress(5, 0)`; that exit is
   `debug_only` (D7). This port used to refuse the FileFilter input with exit 8,
   blaming the zero extent; that refusal is gone, and the test that recorded
   the difference passes as
   `a_short_input_never_reaches_the_seed_loop_as_in_the_cpp_release_build`.
   The derived inputs are built in the tests by the recorded rules and pinned
   to the executed files by SHA-1.
3. **A changed isotope abundance finds the intended features** (`CPP-247`, the
   one designed difference of the algorithm). The tool passes
   `Options::default()` to the algorithm (`src/cli/tools/feature_finder_centroided.rs`),
   whose `AbundanceOverride::Intended` computes the two-isotope override the
   source intends instead of the source's stray-peak override. With
   `-algorithm:isotopic_pattern:abundance_12C 90` on FeatureFinderCentroided_1
   the executed C++ Release tool finds 0 seeds, 0 candidates and 0 features
   (case `ffc1_abundance_12C_90`, three repetitions), and this port finds 18
   seeds, one candidate and one feature, exit 0 — the values of the adapted
   Release replay of the intended override
   (`a_changed_abundance_finds_the_intended_features_where_the_cpp_release_build_finds_none`).
4. **A NaN FAIMS voltage is refused** with exit 6 rather than entering an
   ordered set that cannot hold it: `FaimsHelper::get_compensation_voltages`
   returns an error, which the framework maps to `Invalid parameter: …`. No
   native reader produces such a spectrum; the source would insert the NaN into
   a `std::set<double>` and break its ordering.
5. **The ion-mobility peak type in the message.** The source prints
   `imPeakTypeToString(spec.getIMPeakType())`; its mzML reader sets `IM_PROFILE`
   on every spectrum with an ion-mobility array, except one carrying
   `MS:1003441` (ion mobility centroid frame), which becomes `IM_CENTROIDED`.
   The native spectrum has no ion-mobility peak type and the native reader does
   not keep `MS:1003441`, so the port always prints `im_profile` — the text both
   executed ion-mobility cases printed. A file with `MS:1003441` would print
   `im_centroided` in C++.
6. **`-log`, `writeDebug_` and the timing line are not ported** by the
   framework, so the C++ `TOPP.log` and the debug dump of the `algorithm:`
   parameters at debug level 3 have no counterpart. The tests still run in a
   temporary working directory, so a later `-log` implementation cannot write
   into the repository.
7. **`updateRanges` is not called.** The source updates the experiment's ranges
   after loading (and logs `Update ranges was called but ranges were already
   up-to-date` on most inputs); native ranges are computed on demand.
8. **The overall convex hull is not expanded.** The source's clean-up also calls
   `ft.getConvexHull().expandToBoundingBox()`, which mutates a cache of the
   union of the mass-trace hulls; the native overall hull is computed on demand
   from the mass-trace hulls, so there is nothing to expand. Upstream removed
   that call after the pin for the same reason ("Avoid mutating the cached
   overall feature hull").
9. **An empty mass-trace hull stays empty.** `ConvexHull2D::expandToBoundingBox`
   would replace it with the four corners of an empty bounding box, which are
   `±DBL_MAX` coordinates; the native hull is left alone. The algorithm never
   produces an empty hull.
10. **The debug listing above level 10** lists metadata in key order with Rust's
   shortest round-trip number formatting, where the source lists keys in
   meta-registry order and formats with `StringUtils::toStr`.
11. **`-threads` sizes the seed loop and cannot change a result.** The wrapper
   passes `ToolContext::thread_policy` as
   `feature_finder_picked::algorithm::Options::threads`, which sizes the
   seed-extension loop — the port's form of the source's single `#pragma omp
   parallel for` (`FeatureFinderAlgorithmPicked.cpp:595`). The loop returns its
   results in seed order and every later step is serial, so the determinism
   contract holds strictly: `-threads 0`, `1`, `2`, `4` and `8` write
   byte-identical files, unique ids included. The source's results are
   schedule-independent for the same reason (its `tmp_feature_map` is keyed by
   seed index), but it reaches that by a shared map written from inside the
   parallel region.
12. **Loader strictness.** The tool asks the mzML reader for the source's
    tolerance of dangling `softwareRef` and `defaultDataProcessingRef`
    (decision D10); everything else stays at the reader's strict library
    defaults, so an input the C++ reader repairs silently in another way can
    still be refused here. `docs/MZML_HEADER_SUPPORT.md` lists what the option
    covers.
13. **Non-finite values in the mzML are refused by the reader.** The native
    mzML reader rejects a decoded NaN or infinite m/z or intensity
    (`nonfinite binary value`) and a non-finite `scan start time`, so the tool
    exits 3 before the algorithm runs. The C++ Release tool reads such values
    and hands them to the algorithm, which the library port follows
    ([FEATURE_FINDER_PICKED_SUPPORT](FEATURE_FINDER_PICKED_SUPPORT.md),
    *Non-finite input*). Executed by the adversarial review of this lane
    (`../oracle/ffap-sem-ver1`, driver `write_mod.cpp` writing FFC_1 with one
    value replaced through the Release `MzMLFile`, two repetitions each,
    identical), with the FFC_1 INI and `-test`:

    | Input | C++ Release | This port |
    |---|---|---|
    | `ffc1_mz0_neginf.mzML`: the first m/z `-inf` | exit 8, `FeatureFinder can only operate on spectra that contain peaks with positive m/z values. …` | exit 3, `Unable to read file (parse error on line 0: nonfinite binary value)` |
    | `ffc1_mz0_posinf.mzML`: the first m/z `+inf` (sorted to the end) | exit 8, `the value '12' was used but is not valid; IsotopeDistribution not precalculated. Maximum allowed index is 0` | exit 3, the same reader error |
    | `ffc1_rtlast_posinf.mzML`: the last `scan start time` `inf` | exit 0, `Found 0 seeds for charge 2.`, an empty map | exit 3, `… nonfinite scan start time` |

    The port's exits were run locally on those three files (release build of
    this branch). Given the same values in memory, the library returns what the
    C++ algorithm returns: `mz_neginf_first`, `mz_posinf_last` and
    `rt_posinf_last` of `non_finite_inputs_match_the_linux_release_build`. The
    reader's strictness is the mzML port's (`docs/MZML_HEADER_SUPPORT.md`), not
    this tool's.

14. **The isotope-window ceiling exits 8 where the C++ tool exits 12.** Below
    `vector::max_size()` the source's step 2.5 allocates `56 * count` bytes,
    which fails or not depending on the memory of the process (FFC_1 with its
    last m/z `2e18`: `std::bad_alloc`, reported by TOPPBase's outer handler
    with exit 12). The algorithm refuses every count above
    `Limits::max_isotope_windows` instead (lead decision D6), with its own
    message, which the tool reports as the other algorithm errors, exit 8. The
    debug files are the executed ones (`tool_2e18`). Above `max_size()` the
    tool exits 12 with the source's text, as the executed tool does.

15. **A seed-loop crash exits 8.** With `feature:min_isotope_fit` 0 a seed
    whose best isotope pattern stayed empty makes the algorithm read past an
    empty vector, and the executed tool dies with SIGSEGV (shell status 139,
    case `avg0`). The algorithm refuses there (lead decision D1), and the tool
    reports the refusal as the other algorithm errors, exit 8. The executed
    process had written its console lines up to the FAIMS line; the
    `Found <n> seeds for charge <c>.` lines it had written to `std::cout` were
    still buffered and are lost with the process, while this port prints them
    before the error. In a debug run the port writes `debug/log.txt` up to the
    last byte the executed file buffer had flushed and the feature files the
    executed process had written (the library-level evidence of
    `feature_finder_picked_instrumentation`); the tool reaches that only where
    no seed reaches the fit first (which terminates it as in case b1).

16. **Features of infinite intensity or width cannot be written.** Finite but
    very large retention times make the algorithm's `float` intensity and
    width overflow, and the C++ Release tool writes such features
    (`<intensity>inf</intensity>`, `FWHM` `inf`) and exits 0. The algorithm
    port computes the same features (the source's non-finite values are
    stored since fix round 4), and the tool prints the same console lines,
    but the native featureXML writer refuses a non-finite feature value, so
    the tool exits 3 (`Unable to read file (parse error on line 0: nonfinite
    feature value)`, the writer's error reported through the framework's
    file-error mapping) and writes no output file. The writer's strictness
    belongs to the featureXML port, not to this tool; lead decision D13 of
    wave 5 splits the writer's refusal and the misleading "Unable to read
    file" wording of a write failure off into a separate task, and this
    difference stays as recorded until that task decides. Executed
    (`../oracle/ffap-complete-fix4`, `node/run_tool.sh`, two runs each,
    identical apart from the timing lines), FFC_1's input with every
    `scan start time` value `v` written as `ve36` or `ve39`, FFC_1 INI and
    `-test`:

    | Input | C++ Release | This port |
    |---|---|---|
    | `rt_e36.mzML` | exit 0, 9 features, every intensity `inf`, finite widths | exit 3, the writer error above |
    | `rt_e39.mzML` | exit 0, 9 features, every intensity and `FWHM` `inf` | exit 3, the same |
    | `rt_e39.mzML`, `-algorithm:feature:rt_shape asymmetric` | exit 0, 8 features, every intensity and `FWHM` `inf` | exit 3, the same |

    `infinite_feature_values_are_refused_by_the_featurexml_writer` pins both
    sides (`rt_scaled_release.tsv`).

## Checked boundaries and evidence

- **Bounded work.** The wrapper itself walks the spectra once for the
  ion-mobility check and once for the FAIMS split; both walks are bounded by
  `FaimsHelper::MAX_SPECTRA` respectively by the loader's own limits. The
  per-group loop runs the algorithm once per detected voltage, and there are at
  most as many voltages as spectra; the accumulated feature map is bounded by
  `FeatureMap::MAX_ITEMS`, and the merge by
  `FeatureOverlapFilter::MAX_FEATURES` and
  `FeatureOverlapFilter::MAX_CANDIDATE_VISITS`. The
  annotation steps are bounded by `FeatureMap::MAX_ITEMS` through
  `set_primary_ms_run_path` and `for_each_unique_id`. Loading, seed loading and
  storing inherit the reader and writer ceilings.
- **Atomicity.** `finish_features` takes the map by value and returns it, so a
  failed annotation leaves no half-annotated map; the output file is written
  last, through `FileHandler::store_feature_map`, which replaces the
  destination only after a successful serialisation.
- **No panics on untrusted input.** Every branch returns an exit code or an
  `Error`; the only indexing is `spectra[0]` after the emptiness check.
- **Evidence.** Tier 1 for every asserted exit code, diagnostic and absent
  output (29 executed C++ cases, each run twice and reproduced; 18 more for the
  FAIMS closure, `../oracle/b11-faims`), for the
  `-write_ini` defaults, which are compared with the executed file line by line
  with exact numbers and as a decoded parameter tree, and for the feature
  output, which is compared decoded against the retained expectation and
  measured against the C++ Release build of `bc9cc12`/`174b576` as the table
  above records. Tier 3 for the registration details. The synthetic inputs are
  derived in the tests from the retained
  `FeatureFinderCentroided_1_input.mzML` by the recorded rules and pinned to the
  executed files by digest, so no derived megabyte enters the repository.
- **No C++ oracle.** The cross-voltage merge, because the C++ path never
  reaches it. It is pinned against the derivation in *What the merge is meant
  to do* and against hand-derived numbers, at the tool level
  (`the_faims_merge_joins_the_two_voltages_of_each_analyte`,
  `three_faims_voltages_collapse_a_cluster_of_three_into_one_feature`) and at
  the library level (`tests/feature_overlap_filter.rs`, the corrected-merge
  section), always beside the executed `c2_*` cases that hold what the source
  does instead. Everything the merge is given — the split, the per-group runs,
  the group order and the annotation — is pinned against executed C++.
- **Not asserted.** The Debug exit codes of the two `debug_only` cases (D7); the
  C++ Release behaviour is asserted instead for both, `FFC_im_arrays_ms2_only`
  and `FFC_FileFilter_44_force`. The `1e-9` comparison
  against the C1 oracle outputs of `FFC_seeds`, `FFC_asymmetric` and
  `FFC_debug5`, and the `-algorithm:fit:max_iterations` boundary sweep, are
  package B10's and are not asserted here.
