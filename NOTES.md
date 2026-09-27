# B6 lane notes (SHARED-01: DefaultParamHandler trait)

Branch: worktree-wf_ee239125-5be-2 (merged prior incarnation worktree-wf_bcdd783d-f3f-2, notes only).

## What is missing
- No trait over the defaults()/from_param()/from_param_with_warnings()/to_param() pattern.
  Pattern source: src/processing/peak_picking.rs ~1853-1985 (PeakPickerHiRes), uses
  crate::param::DefaultParamHandler (src/param/handler.rs) set_defaults + set_parameters_with.
  PEAK_PICKER_HI_RES_NAME const is the handler name.
- Ledger targets (docs/core-sdk-coverage.json, scope says param gap):
  PROCESSING partial: PeakPickerIterative (src/processing/iterative.rs), NLargest, ThresholdMower,
  Normalizer, LinearResamplerAlign (all src/processing.rs), SignalToNoiseEstimatorMeanIterative
  (src/processing/mean_noise.rs), GaussFilter, SavitzkyGolayFilter (src/processing/smoothing.rs);
  complete (C++ .cpp names DefaultParamHandler only): WindowMower, RankScaler, SqrtScaler.
  CHEMISTRY partial: TheoreticalSpectrumGenerator (src/chemistry/theoretical.rs), SpectrumAnnotator
  (src/chemistry/spectrum_annotator.rs), NucleicAcidSpectrumGenerator.
  Also mentions: ANALYSIS/OPENSWATH/PeakPickerChromatogram (not mine).
- SignalToNoiseEstimatorMedian already has defaults/from_param in src/processing/peak_picking/noise.rs.

## Evidence
- C++ getDefaults: git -C ../OpenMS4-tests/packages/core show bc9cc12:src/openms/source/<...>.cpp
- class tests: src/tests/class_tests/openms/source/<Name>_test.cpp
- Oracle: ../oracle/B6-param-defaults/ on ibminode06, Release install, dump getDefaults via ParamXMLFile.

## Done
- Oracle ../oracle/B6-param-defaults RAN OK on ibminode06 Release install (results/, no DIFFER).
  attempt1/ = first try, failed: SignalToNoiseEstimatorMeanIterative<MSChromatogram> does not compile
  (getMZ on ChromatogramPeak, header :262,:279) -> C++ defect candidate; case dropped.
  Unknown-param warning in C++: "Warning: X received the unknown parameter 'k'!" (Rust handler says
  "X: unknown parameter 'k'" -- left as is, FeatureFinderPicked maps it itself).
- Trait `DefaultParamHandling` in src/param/handler.rs (re-exported from param): NAME, SUBSECTIONS,
  CHECK_DEFAULTS, WARN_EMPTY_DEFAULTS, defaults(), to_param(), update_members(), keep_native_members(),
  provided handler/from_param/_with_warnings/_with_handler/set_parameters.
  Adopted on PeakPickerHiRes and SignalToNoiseEstimatorMedian.
- C++ sources readable at OpenMS4-R/.reference/openms4-core-bc9cc12 (main checkout). The worktree
  sandbox refuses a git invocation against another repo; use Read/cat with literal paths.

- 11 PROCESSING impls done (src/processing/parameters.rs, window_mower.rs, smoothing.rs,
  mean_noise.rs, iterative.rs); tests/default_param_handling.rs 12/12 pass on spock (tier 1 vs
  tests/data/param_defaults/*.ini copied from oracle results; hashes in scratchpad B6/fixture_hashes.txt).
  Param::replace_value added (src/param.rs); helpers unsigned_parameter/int_parameter/flag_parameter/
  set_flag in param::handler. GaussFilter gained inactive_width + write_log_messages; MeanIterative
  gained range_parameters; SavitzkyGolay keeps frame_parameter.

## Next step
1. Chemistry: TheoreticalSpectrumGenerator (Rust intensities are f32! source double -> to_param
   precision problem; decide), SpectrumAnnotator, NucleicAcidSpectrumGenerator.
2. Docs: docs/DEFAULT_PARAM_HANDLING_SUPPORT.md, provenance json, doc coverage 100% for touched
   modules (param.rs 22%, handler.rs 54%, mean_noise now better), CI slice, ledger, report.
