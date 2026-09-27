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

## Next step
1. Design trait in src/param/handler.rs (DefaultParameters? name TBD) -- write it, commit.
2. Implement for PeakPickerHiRes (reference), then the 14 headers.
3. Oracle for default Param dumps.
