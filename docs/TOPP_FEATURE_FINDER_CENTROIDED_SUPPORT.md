# FeatureFinderCentroided: supported source subset

The TOPP tool that detects two-dimensional features in centroided LC-MS data.
Ported by package **C5-FFC-WRAPPER** of the early TOPP bundle from
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
| FAIMS input | **refused** with exit 11 (decision D5); the closure is package B11's |
| The picked algorithm | ported up to seed selection; every accepted input ends with `Error: unsupported: …` and exit 11 until package B7 lands |
| Primary MS run path, unique ids, `QUANTITATION` processing record, hull and subordinate clean-up, featureXML store | implemented; tested directly, because no run reaches them yet |
| `TOPP_FeatureFinderCentroided_1` output, `-seeds` results, `feature:rt_shape asymmetric`, `-debug 5`, thread counts | package B10 |

The tool therefore stays **partial** in the ledger: it decides every input the
C++ tool decides before the algorithm, and produces no feature map.

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
| `registerStringOption_("faims_merge_features", "<true/false>", "true", …, false)` + `setValidStrings_` | `ToolSpec::register_string_option` + `set_valid_strings`; registered and validated, without effect while FAIMS input is refused |
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
| `IMDataConverter::splitByFAIMSCV(std::move(exp))`, `has_faims`, the per-group loop, the seed filter, the FAIMS_CV annotation and `FeatureOverlapFilter::mergeFAIMSFeatures` | not ported here (decision D5): `FaimsHelper::get_compensation_voltages` decides whether the input is FAIMS input, and a non-empty set is refused with `FeatureFinderCentroided::faims_refusal_message`. `crate::kernel::im_data_converter` (B8) and `crate::processing::feature_overlap_filter` (B9) hold the library halves; package B11 joins them |
| `OPENMS_LOG_INFO << "Not FAIMS compensation voltages found …"` (in `splitByFAIMSCV`) | `FeatureFinderCentroided::NO_FAIMS_MESSAGE` on the output stream |
| `OPENMS_LOG_INFO << "FAIMS data detected with N compensation voltage(s)."` | `FeatureFinderCentroided::faims_detected_message` |
| `OPENMS_LOG_INFO << "Processing FAIMS CV group: …"`, `"Combined N features …"`, `"FAIMS feature merge: …"` | not reached: the refusal comes first |
| `FeatureFinderAlgorithmPicked ff; ff.run(std::move(group), features_cv, feafi_param, seeds_cv)` | `crate::analysis::feature_finder_picked::algorithm::run_with_options`, called once, because a refused FAIMS input is the only case with more than one group |
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
case (C1 = `../oracle/topp-early-bundle`, C5 = `../oracle/ffc-wrapper-c5`).

| Input | Exit | Diagnostic | Oracle case |
|---|---|---|---|
| Registered workflow `TOPP_FeatureFinderCentroided_1` | C++ 0; this port 11 until B7 | `Error: unsupported: FeatureFinderAlgorithmPicked seed extension, …` | C1 `TOPP_FeatureFinderCentroided_1` |
| `-write_ini` (with and without `-test`) | 0 | — | C1 `FFC_write_ini`, `FFC_write_ini_test` |
| MS2 spectra only | 4 | `Error: File empty (the file 'Error: No MS1 spectra in input file.' is empty)` | C1 `FFC_ms2_only` |
| Per-peak ion mobility, with or without a unit | 11 | `Error: Input contains per-peak ion mobility data (IM_PEAK, im_profile) …` | C1 `FFC_im_peak_with_units`, `FFC_im_peak_without_units` |
| Per-peak ion mobility **and** profile data | 11, the ion-mobility message | the check order of `main_` | C5 `c5_im_peak_profile_noforce` |
| Ion-mobility arrays on MS2 spectra only | not refused | — | C1 `FFC_im_arrays_ms2_only` (the C++ exit is Debug-only) |
| First spectrum stored profile, no `-force` | 8 | `Error: Unexpected internal error (Error: Profile data provided but centroided spectra expected. …)` | C1 `FFC_profile_noforce`, `FFC_FileFilter_44_noforce`, C5 `c5_first_profile_only` |
| First spectrum profile, `-force` | not refused | — | C1 `FFC_profile_force` |
| Profile term before `MS:1000525`, no `-force` | not refused: the reader resets the type | — | C1 `FFC_profile_then_spectrum_representation` |
| Every spectrum but the first stored profile | not refused: only `exp[0]` is checked | — | C5 `c5_later_profile_only` |
| Every MS1 peak below the intensity range | 8 | `Error: Unexpected internal error (FeatureFinder needs updated ranges on input map. Aborting.)` | C5 `c5_negative_intensities` |
| `-seeds` that is not featureXML | 6 | `Input file '…' has invalid format 'mzML'. Valid formats are: 'featureXML'.` | C5 `c5_seeds_not_featurexml` |
| FAIMS input with an unreadable `-seeds` file | 3 | `Error: Unable to read file (…)`, before any FAIMS message | C5 `c5_faims_corrupt_seeds` |
| Any FAIMS input (one voltage, two voltages, a voltage on half the spectra, `-faims_merge_features false`) | C++ 8; this port 11 | `Error: FAIMS input is not supported by this port of FeatureFinderCentroided yet (compensation voltages … V): …` | C1 `FFC_faims_*`, C5 `c5_faims_partial_cv` |
| FAIMS **profile** input without `-force` | 8, the profile message | the profile check precedes the split | C1 `FFC_faims_interleaved_noforce` |
| `-algorithm:feature:rt_shape bogus` | 6 | `Invalid string parameter value 'bogus' … Valid values are: 'symmetric,asymmetric'.` | C1 `FFC_invalid_rt_shape` |
| `-out` without an extension | accepted | — | C1 `FFC_out_no_extension` |

No branch writes `-out` before the store, and every refusal above was checked
to leave no output file.

## Preserved source conventions

- **Check order.** Load, empty input, per-peak ion mobility, profile, seeds,
  `algorithm:` parameters, FAIMS, algorithm, annotation, clean-up, store. Three
  of these orders are observable and each is pinned by an executed case:
  ion mobility before profile, profile before FAIMS, seeds before FAIMS.
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

1. **FAIMS input is refused** (decision D5). The C++ tool splits by compensation
   voltage and fails afterwards: its groups are built with `addSpectrum`, which
   leaves them without per-MS-level ranges, so `FeatureFinderAlgorithmPicked`
   throws `the value '1' was used but is not valid; No ranges for this MS level`
   and every FAIMS input exits 8. Even with that fixed, `mergeFAIMSFeatures`
   removes every feature, because the features still carry unique id 0 when it
   keys removal by id. Reproducing either would be pointless, and pooling the
   voltages silently would be wrong, so the port refuses such input with exit 11
   and writes nothing. Package B11 ports the closure and removes the refusal.
2. **The algorithm stops after seed selection** (package B7). Every input that
   passes the wrapper ends with exit 11 and the message quoted above. Two
   further observed endings belong to the algorithm, not the wrapper, and are
   recorded here only so a reader is not surprised: the numpress profile fixture
   with `-force` ends with `Error: Unexpected internal error
   (FeatureFinderAlgorithmPicked needs a retention-time and an m/z range of
   positive width …)` and exit 8, where the C++ Debug build fails an intensity
   precondition; `FileConverter_31_output.mzML` ends with the not-ported message
   and exit 11, where the C++ Debug build fails a `ProgressLogger` precondition.
   Both C++ exits are Debug-only and are not expectations.
3. **Two `algorithm:` values the ported stage refuses** reach the user through
   this tool: `-algorithm:write_debug true`, whose source debug output reads an
   undeclared parameter and writes into the working directory, and a non-default
   `-algorithm:isotopic_pattern:abundance_12C` or `abundance_14N`, whose source
   handling keeps a stray peak in the isotope distribution. Both end in exit 11
   with the algorithm's message; `FEATURE_FINDER_PICKED_SUPPORT.md` records why.
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
11. **`-threads` changes nothing yet.** The source parallelises only seed
   extension (`FeatureFinderAlgorithmPicked.cpp:595`), which is not ported; the
   ported stage is serial in the source too. The entry point takes no thread
   policy, so the wrapper passes none; B7 adds the parameter and B10 passes
   `ToolContext::thread_policy`. The determinism contract is met trivially: no
   result depends on a thread count.
12. **Loader strictness.** The tool asks the mzML reader for the source's
    tolerance of dangling `softwareRef` and `defaultDataProcessingRef`
    (decision D10); everything else stays at the reader's strict library
    defaults, so an input the C++ reader repairs silently in another way can
    still be refused here. `docs/MZML_HEADER_SUPPORT.md` lists what the option
    covers.

## Checked boundaries and evidence

- **Bounded work.** The wrapper itself walks the spectra once for the
  ion-mobility check and once for the FAIMS voltages; both walks are bounded by
  `FaimsHelper::MAX_SPECTRA` respectively by the loader's own limits. The
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
  output (29 executed C++ cases, each run twice and reproduced), and for the
  `-write_ini` defaults, which are compared with the executed file line by line
  with exact numbers and as a decoded parameter tree. Tier 3 for the
  registration details and the annotation steps, whose end-to-end comparison is
  package B10's. The synthetic inputs are derived in the tests from the retained
  `FeatureFinderCentroided_1_input.mzML` by the recorded rules and pinned to the
  executed files by digest, so no derived megabyte enters the repository.
- **Not asserted.** The exit codes of the two cases whose C++ runs end in a
  Debug-only precondition, and anything about the feature map a complete run
  would produce.
